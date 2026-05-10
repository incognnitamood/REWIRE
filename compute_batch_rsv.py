"""
Batch RSV computation for all completed drug graphs.

Run with:
    python compute_batch_rsv.py
"""

from __future__ import annotations

import sys
import pickle
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from tqdm import tqdm

# Add repo src paths for imports
REPO_ROOT = Path(__file__).parent / "REWIRE"
SRC_ROOT = REPO_ROOT / "src"
SRC_P2 = SRC_ROOT / "P2"

sys.path.insert(0, str(SRC_ROOT))
sys.path.insert(0, str(SRC_P2))

from ppi_graph import build_graph
from rsv_compute import compute_rsv


def _load_drug_ids(path: Path) -> List[str]:
    if not path.exists():
        raise FileNotFoundError(f"drug_ids_completed.txt not found: {path}")
    with path.open("r", encoding="utf-8") as f:
        return [line.strip() for line in f if line.strip()]


def _normalize_key(name: str) -> str:
    return name.strip().lower()


def _load_manifest() -> Tuple[Dict[str, Dict[str, object]], str]:
    """Load a manifest mapping from drug_id/drug_name to targets and display name.

    Expected primary location: src/P2/manifest.csv
    Fallbacks: src/P2/drug_graphs/manifest.csv, src/P2/drug_50.csv
    """

    candidates = [
        SRC_P2 / "manifest.csv",
        SRC_P2 / "drug_graphs" / "manifest.csv",
        SRC_P2 / "drug_50.csv",
    ]

    manifest_path = None
    df = None

    for path in candidates:
        if not path.exists():
            continue
        try:
            probe = pd.read_csv(path)
        except Exception:
            continue

        probe.columns = [c.strip().lower() for c in probe.columns]
        has_target_proteins = "target_proteins" in probe.columns or "targets" in probe.columns
        has_target_gene = "target_gene" in probe.columns or "gene_symbol" in probe.columns

        # Skip manifests that don't provide targets.
        if not (has_target_proteins or has_target_gene):
            continue

        df = probe
        manifest_path = path
        break

    if df is None or manifest_path is None:
        raise FileNotFoundError(
            "No manifest CSV with targets found in expected locations."
        )

    has_target_proteins = "target_proteins" in df.columns or "targets" in df.columns
    has_target_gene = "target_gene" in df.columns or "gene_symbol" in df.columns
    has_drug_id = "drug_id" in df.columns
    has_drug_name = "drug_name" in df.columns

    if not (has_drug_id or has_drug_name):
        raise ValueError(f"Manifest missing drug_id or drug_name columns: {manifest_path}")

    manifest_map: Dict[str, Dict[str, object]] = {}

    if has_target_proteins:
        target_col = "target_proteins" if "target_proteins" in df.columns else "targets"
        id_col = "drug_id" if has_drug_id else "drug_name"

        for _, row in df.iterrows():
            drug_key = str(row[id_col]).strip()
            drug_name = str(row["drug_name"]).strip() if has_drug_name else drug_key
            targets_raw = str(row.get(target_col, "")).strip()
            targets = [t.strip() for t in targets_raw.split(",") if t.strip()]

            manifest_map[_normalize_key(drug_key)] = {
                "drug_id": drug_key,
                "drug_name": drug_name,
                "targets": targets,
                "pkl_file": row.get("pkl_file", None),
            }
    elif has_target_gene:
        # Long format: one row per drug-target
        id_col = "drug_id" if has_drug_id else "drug_name"
        target_col = "target_gene" if "target_gene" in df.columns else "gene_symbol"

        grouped = df.groupby(id_col)[target_col].apply(list).to_dict()
        for drug_key, targets in grouped.items():
            drug_key_str = str(drug_key).strip()
            # Try to find a display name if available
            drug_name = drug_key_str
            if has_drug_name:
                try:
                    drug_name = df.loc[df[id_col] == drug_key, "drug_name"].iloc[0]
                except Exception:
                    pass

            manifest_map[_normalize_key(drug_key_str)] = {
                "drug_id": drug_key_str,
                "drug_name": str(drug_name).strip(),
                "targets": [str(t).strip() for t in targets if str(t).strip()],
                "pkl_file": None,
            }
    else:
        raise ValueError(
            f"Manifest missing target_proteins/targets or target_gene columns: {manifest_path}"
        )

    return manifest_map, str(manifest_path)


def _resolve_pkl_path(drug_id: str, manifest_entry: Dict[str, object] | None) -> Path:
    base_dir = SRC_P2 / "drug_graphs"
    candidates = []

    if manifest_entry and manifest_entry.get("pkl_file"):
        candidates.append(str(manifest_entry["pkl_file"]))

    candidates.extend([
        f"{drug_id}.pkl",
        f"{drug_id.lower()}.pkl",
        f"{drug_id.lower().replace(' ', '_')}.pkl",
    ])

    for name in candidates:
        path = base_dir / name
        if path.exists():
            return path

    return base_dir / f"{drug_id}.pkl"


def _filter_targets(targets: List[str], G0) -> List[str]:
    if not targets:
        return []

    missing = [t for t in targets if t not in G0]
    if missing:
        print(f"[WARN] Missing targets not in G0: {missing[:5]}" + (" ..." if len(missing) > 5 else ""))

    return [t for t in targets if t in G0]


def main() -> int:
    ppi_path = SRC_P2 / "ppi_genes.csv"
    drug_ids_path = SRC_P2 / "drug_graphs" / "drug_ids_completed.txt"

    print(f"Loading G0 from {ppi_path}")
    G0 = build_graph(str(ppi_path))
    if G0 is None:
        print("ERROR: Failed to load baseline PPI graph.")
        return 1

    drug_ids = _load_drug_ids(drug_ids_path)
    if not drug_ids:
        print("ERROR: No drug IDs found.")
        return 1

    manifest_map, manifest_path = _load_manifest()
    print(f"Loaded manifest: {manifest_path}")

    rsv_vectors: List[List[float]] = []
    out_drug_ids: List[str] = []
    out_drug_names: List[str] = []
    failed_drugs: List[str] = []

    partial_path = Path(__file__).parent / "drug_rsv_partial.npy"

    for idx, drug_id in enumerate(tqdm(drug_ids, desc="Computing RSV", unit="drug"), start=1):
        try:
            entry = manifest_map.get(_normalize_key(drug_id))
            if entry is None:
                print(f"[WARN] No manifest entry for drug_id={drug_id}. Using empty targets.")
                entry = {"drug_id": drug_id, "drug_name": drug_id, "targets": [], "pkl_file": None}

            pkl_path = _resolve_pkl_path(drug_id, entry)
            if not pkl_path.exists():
                raise FileNotFoundError(f"PKL not found for {drug_id}: {pkl_path}")

            with pkl_path.open("rb") as f:
                G_drug = pickle.load(f)

            targets = entry.get("targets", [])
            if isinstance(targets, str):
                targets = [t.strip() for t in targets.split(",") if t.strip()]

            targets = _filter_targets(list(targets), G0)

            rsv = compute_rsv(G0, G_drug, targets)
            vector = rsv.get("vector") if isinstance(rsv, dict) else rsv

            if not vector or len(vector) != 4:
                raise ValueError(f"RSV vector invalid for {drug_id}: {vector}")

            rsv_vectors.append([float(v) for v in vector])
            out_drug_ids.append(str(entry.get("drug_id", drug_id)))
            out_drug_names.append(str(entry.get("drug_name", drug_id)))

        except Exception as exc:
            print(f"[WARN] Failed drug {drug_id}: {exc}")
            failed_drugs.append(drug_id)

        if idx % 10 == 0:
            if rsv_vectors:
                np.save(partial_path, np.asarray(rsv_vectors, dtype=np.float32))
                print(f"[INFO] Saved partial matrix to {partial_path} (rows={len(rsv_vectors)})")
            else:
                np.save(partial_path, np.zeros((0, 4), dtype=np.float32))

    # Save final outputs
    output_dir = Path(__file__).parent
    matrix_path = output_dir / "drug_rsv_matrix.npy"
    ids_path = output_dir / "drug_ids.txt"
    names_path = output_dir / "drug_names.txt"

    matrix = np.asarray(rsv_vectors, dtype=np.float32)
    np.save(matrix_path, matrix)

    with ids_path.open("w", encoding="utf-8") as f:
        for item in out_drug_ids:
            f.write(f"{item}\n")

    with names_path.open("w", encoding="utf-8") as f:
        for item in out_drug_names:
            f.write(f"{item}\n")

    print("\nSuccessfully processed: {}/{} drugs".format(len(out_drug_ids), len(drug_ids)))
    print(f"Failed drugs: {failed_drugs}")
    print(f"Matrix shape: {matrix.shape}")
    print(f"Saved drug_rsv_matrix.npy")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
