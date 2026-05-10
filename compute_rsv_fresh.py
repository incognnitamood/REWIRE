"""Compute RSVs by simulating binding fresh (no PKL reuse).

Run with:
    python compute_rsv_fresh.py
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
from rsv_compute import compute_rsv


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


def _first_target_with_neighbor(G0, targets: List[str]) -> str | None:
    for t in targets:
        if t in G0 and len(list(G0.neighbors(t))) > 0:
            return t
    return None


def _target_subgraph(G: object, targets: List[str]) -> object:
    nodes = set()
    for t in targets:
        if t not in G:
            continue
        nodes.add(t)
        nodes.update(G.neighbors(t))
    return G.subgraph(nodes).copy() if nodes else G


def main() -> int:
    ppi_path = SRC_P2 / "ppi_genes.csv"
    drug_csv = SRC_P2 / "drug_50.csv"

    print(f"Loading G0 from {ppi_path}")
    G0 = build_graph(str(ppi_path))
    if G0 is None:
        print("ERROR: Failed to load baseline PPI graph.")
        return 1

    if not drug_csv.exists():
        print(f"ERROR: drug_50.csv not found at {drug_csv}")
        return 1

    df = pd.read_csv(drug_csv)
    df.columns = [c.strip().lower() for c in df.columns]

    # Accept the requested schema or fall back to minimal columns.
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

    rsv_vectors = []
    drug_ids = []
    drug_names = []

    warning_count = 0
    any_weight_changed = False

    for _, row in tqdm(df.iterrows(), total=len(df), desc="Computing RSV", unit="drug"):
        drug_id = str(row.get("drug_id", "")).strip()
        drug_name = str(row.get("drug_name", drug_id)).strip()

        targets = _parse_list(str(row.get("target_proteins", "")))
        affinities = _parse_floats(str(row.get("binding_affinities", "")))

        if not targets:
            print(f"WARNING: drug {drug_name} has no targets. Skipping.")
            warning_count += 1
            continue

        if affinities and len(affinities) != len(targets):
            print(f"WARNING: drug {drug_name} has mismatched targets/affinities. Using defaults.")
            affinities = []

        affinity_map = {t: a for t, a in zip(targets, affinities)} if affinities else {}

        try:
            G_drug, _failed = simulate_binding(G0, targets, affinity_map)
        except Exception as exc:
            print(f"WARNING: drug {drug_name} — simulate_binding failed: {exc}")
            warning_count += 1
            continue

        # Verify perturbation on the first valid target
        probe_target = _first_target_with_neighbor(G0, targets)
        if probe_target is None:
            print(f"WARNING: drug {drug_name} — no valid targets found in G0")
            warning_count += 1
        else:
            neighbor = next(iter(G0.neighbors(probe_target)))
            w0 = G0[probe_target][neighbor].get("weight")
            wd = G_drug[probe_target][neighbor].get("weight")

            if w0 == wd:
                print(f"WARNING: drug {drug_name} — simulation had no effect")
                warning_count += 1
            else:
                print(f"OK: drug {drug_name} — weight changed from {w0} to {wd}")
                any_weight_changed = True

        try:
            G0_sub = _target_subgraph(G0, targets)
            Gd_sub = _target_subgraph(G_drug, targets)
            rsv = compute_rsv(G0_sub, Gd_sub, targets)
            vector = rsv.get("vector") if isinstance(rsv, dict) else rsv
            if not vector or len(vector) != 4:
                print(f"WARNING: drug {drug_name} — invalid RSV vector")
                warning_count += 1
                continue

            rsv_vectors.append([float(v) for v in vector])
            drug_ids.append(drug_id)
            drug_names.append(drug_name)
        except Exception as exc:
            print(f"WARNING: drug {drug_name} — compute_rsv failed: {exc}")
            warning_count += 1
            continue

    if not any_weight_changed:
        print(
            "CRITICAL: simulate_binding is not modifying any weights.\n"
            "Check that ppi_graph.py simulate_binding function actually\n"
            "modifies edges and does not return an unmodified copy."
        )
        return 1

    output_dir = Path(__file__).parent
    matrix_path = output_dir / "drug_rsv_matrix.npy"
    ids_path = output_dir / "drug_ids.txt"
    names_path = output_dir / "drug_names.txt"

    matrix = np.asarray(rsv_vectors, dtype=np.float32)
    np.save(matrix_path, matrix)

    with ids_path.open("w", encoding="utf-8") as f:
        for item in drug_ids:
            f.write(f"{item}\n")

    with names_path.open("w", encoding="utf-8") as f:
        for item in drug_names:
            f.write(f"{item}\n")

    print(f"\nSuccessful: {len(drug_ids)} / {len(df)}")
    print(f"Warnings (no perturbation): {warning_count} drugs")
    print(f"Matrix shape: {matrix.shape}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
