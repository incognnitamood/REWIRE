"""Rank drugs by cosine similarity to disease DRS vectors.

Run with:
    python rank_drugs.py
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Dict

import numpy as np


REPO_ROOT = Path(__file__).parent

DRUG_RSV_PATH = REPO_ROOT / "drug_rsv_matrix.npy"
DRUG_NAMES_PATH = REPO_ROOT / "drug_names.txt"
DISEASE_DRS_PATH = REPO_ROOT / "disease_drs.npy"
DISEASE_NAMES_PATH = REPO_ROOT / "disease_names.txt"


def _load_names(path: Path) -> List[str]:
    if not path.exists():
        raise FileNotFoundError(f"Names file not found: {path}")
    with path.open("r", encoding="utf-8") as f:
        return [line.strip() for line in f if line.strip()]


def _find_disease_index(disease_name: str, disease_names: List[str]) -> int | None:
    query = disease_name.strip().lower()
    if not query:
        return None

    # Exact match (case-insensitive)
    for i, name in enumerate(disease_names):
        if name.lower() == query:
            return i

    # Partial match (case-insensitive)
    for i, name in enumerate(disease_names):
        if query in name.lower():
            return i

    return None


def _normalize_matrix(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mean = matrix.mean(axis=0)
    std = matrix.std(axis=0)
    std[std < 1e-8] = 1.0
    return (matrix - mean) / std, mean, std


def rank_for_disease(disease_name: str, top_k: int = 10) -> List[Dict[str, object]]:
    idx = _find_disease_index(disease_name, DISEASE_NAMES)
    if idx is None:
        print("Disease not found. Available diseases:")
        for name in DISEASE_NAMES:
            print(f"- {name}")
        return []

    # Compress ΔEN variance — log1p pulls extreme values (0 or 1) toward center
    rsv_matrix = DRUG_RSV.copy().astype(float)
    rsv_matrix[:, 3] = np.log1p(rsv_matrix[:, 3])
    drs_vector = DISEASE_DRS[idx].copy().astype(float)
    drs_vector[3] = np.log1p(drs_vector[3])

    rsv_norm, rsv_mean, rsv_std = _normalize_matrix(rsv_matrix)
    drs_norm = (drs_vector - rsv_mean) / rsv_std

    weights = np.array([0.15, 0.40, 0.25, 0.20])
    rsv_weighted = rsv_norm * weights
    drs_weighted = drs_norm * weights

    drs_norm_val = np.linalg.norm(drs_weighted) + 1e-8
    scores = []
    for i in range(rsv_weighted.shape[0]):
        rsv_norm_val = np.linalg.norm(rsv_weighted[i]) + 1e-8
        score = float(np.dot(rsv_weighted[i], drs_weighted) / (rsv_norm_val * drs_norm_val))
        scores.append(score)
    scores = np.array(scores)

    order = np.argsort(scores)[::-1]
    top_idx = order[:top_k]

    results = []
    for rank, i in enumerate(top_idx, start=1):
        results.append({
            "rank": rank,
            "drug_name": DRUG_NAMES[i],
            "score": float(scores[i])
        })

    return results


def _print_table(disease_name: str, results: List[Dict[str, object]]) -> None:
    print(f"Top {len(results)} candidates for {disease_name}")
    print("=" * 41)
    print(f"{'Rank':<5} {'Drug Name':<20} {'Score'}")

    for row in results:
        print(f"{row['rank']:<5} {row['drug_name']:<20} {row['score']:.4f}")
    print("")


if __name__ == "__main__":
    DRUG_RSV = np.load(DRUG_RSV_PATH)
    DRUG_NAMES = _load_names(DRUG_NAMES_PATH)
    DISEASE_DRS = np.load(DISEASE_DRS_PATH)
    DISEASE_NAMES = _load_names(DISEASE_NAMES_PATH)

    # Parkinson's Disease demo
    results = rank_for_disease("Parkinson's Disease", top_k=10)
    if results:
        _print_table("Parkinson's Disease", results)

    # All diseases
    for disease in DISEASE_NAMES:
        results = rank_for_disease(disease, top_k=10)
        if results:
            _print_table(disease, results)
