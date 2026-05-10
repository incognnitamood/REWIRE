"""Patch columns 2 and 3 in drug_rsv_matrix.npy by recomputing them.

Run with:
    python patch_rsv_columns.py
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd
from tqdm import tqdm

# Add repo src paths for imports
REPO_ROOT = Path(__file__).parent / "REWIRE"
SRC_ROOT = REPO_ROOT / "src"
SRC_P2 = SRC_ROOT / "P2"

sys.path.insert(0, str(SRC_ROOT))
sys.path.insert(0, str(SRC_P2))

from ppi_graph import build_graph, simulate_binding
from rsv_compute import compute_spectral_gap, compute_entropy_delta


def _parse_list(value: str) -> List[str]:
    if not isinstance(value, str):
        return []
    return [v.strip() for v in value.split(",") if v.strip()]


def _parse_floats(value: str) -> List[float]:
    if not isinstance(value, str):
        return []
    result = []
    for v in value.split(","):
        v = v.strip()
        if not v:
            continue
        try:
            result.append(float(v))
        except ValueError:
            continue
    return result


def main() -> int:
    matrix_path = Path(__file__).parent / "drug_rsv_matrix.npy"
    ids_path = Path(__file__).parent / "drug_ids.txt"
    names_path = Path(__file__).parent / "drug_names.txt"
    drug_csv = SRC_P2 / "drug_50.csv"
    ppi_path = SRC_P2 / "ppi_genes.csv"

    if not matrix_path.exists():
        print(f"ERROR: Missing {matrix_path}")
        return 1
    if not ids_path.exists() or not names_path.exists():
        print("ERROR: Missing drug_ids.txt or drug_names.txt")
        return 1
    if not drug_csv.exists():
        print(f"ERROR: Missing drug_50.csv at {drug_csv}")
        return 1

    matrix = np.load(matrix_path)
    drug_ids = [line.strip() for line in ids_path.read_text(encoding="utf-8").splitlines() if line.strip()]
    drug_names = [line.strip() for line in names_path.read_text(encoding="utf-8").splitlines() if line.strip()]

    if matrix.shape[0] != len(drug_ids):
        print("ERROR: Matrix row count does not match drug_ids.txt")
        return 1

    print(f"Loading G0 from {ppi_path}")
    G0 = build_graph(str(ppi_path))
    if G0 is None:
        print("ERROR: Failed to load baseline PPI graph.")
        return 1

    df = pd.read_csv(drug_csv)
    df.columns = [c.strip().lower() for c in df.columns]

    required_cols = {"drug_id", "drug_name", "target_proteins", "binding_affinities"}
    has_full = required_cols.issubset(set(df.columns))

    if not has_full:
        if {"drug_name", "target_gene"}.issubset(set(df.columns)):
            df = df.rename(columns={"target_gene": "target_proteins"})
            df["binding_affinities"] = ""
            df["drug_id"] = df["drug_name"]
        else:
            print("ERROR: drug_50.csv missing expected columns.")
            print("Expected either: drug_id, drug_name, target_proteins, binding_affinities")
            print("or: drug_name, target_gene")
            return 1

    # Map drug_id -> row index
    id_to_index = {drug_id: i for i, drug_id in enumerate(drug_ids)}

    # Print before values for Imatinib if present
    if "Imatinib" in drug_names:
        idx = drug_names.index("Imatinib")
        print(f"Imatinib BEFORE col2,col3: {matrix[idx, 2]}, {matrix[idx, 3]}")

    for _, row in tqdm(df.iterrows(), total=len(df), desc="Patching RSV", unit="drug"):
        drug_id = str(row.get("drug_id", "")).strip()
        if drug_id not in id_to_index:
            continue

        targets = _parse_list(str(row.get("target_proteins", "")))
        affinities = _parse_floats(str(row.get("binding_affinities", "")))
        if affinities and len(affinities) != len(targets):
            affinities = []

        affinity_map = {t: a for t, a in zip(targets, affinities)} if affinities else {}

        try:
            G_drug, _failed = simulate_binding(G0, targets, affinity_map)
            matrix[id_to_index[drug_id], 2] = compute_spectral_gap(G0, G_drug)
            matrix[id_to_index[drug_id], 3] = compute_entropy_delta(G0, G_drug, targets)
        except Exception as exc:
            print(f"WARNING: {drug_id} failed: {exc}")
            continue

    np.save(matrix_path, matrix)

    if "Imatinib" in drug_names:
        idx = drug_names.index("Imatinib")
        print(f"Imatinib AFTER col2,col3: {matrix[idx, 2]}, {matrix[idx, 3]}")

    print("Patched drug_rsv_matrix.npy")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
