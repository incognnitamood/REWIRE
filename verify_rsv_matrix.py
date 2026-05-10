"""Verify the RSV matrix and report summary statistics.

Run with:
    python verify_rsv_matrix.py
"""

from __future__ import annotations

from pathlib import Path
import numpy as np


def _load_names(path: Path) -> list[str]:
    if not path.exists():
        raise FileNotFoundError(f"drug_names.txt not found: {path}")
    with path.open("r", encoding="utf-8") as f:
        return [line.strip() for line in f if line.strip()]


def _print_col_stats(matrix: np.ndarray, labels: list[str]) -> None:
    print("\nColumn stats:")
    for i, label in enumerate(labels):
        col = matrix[:, i]
        print(f"  {label}:")
        print(f"    min={np.nanmin(col):.6f}  max={np.nanmax(col):.6f}")
        print(f"    mean={np.nanmean(col):.6f} std={np.nanstd(col):.6f}")


def _print_specific_drugs(matrix: np.ndarray, names: list[str]) -> None:
    targets = ["Imatinib", "Nilotinib", "Metformin", "Bexarotene"]
    name_to_index = {name: i for i, name in enumerate(names)}

    print("\nSelected drug RSVs:")
    for name in targets:
        if name not in name_to_index:
            print(f"  {name} RSV: not found")
            continue
        rsv = matrix[name_to_index[name]]
        print(f"  {name} RSV: [{rsv[0]:.3f}, {rsv[1]:.3f}, {rsv[2]:.3f}, {rsv[3]:.3f}]")


def _scan_warnings(matrix: np.ndarray, names: list[str]) -> bool:
    has_warning = False

    for i, row in enumerate(matrix):
        drug_name = names[i] if i < len(names) else f"row_{i}"

        if np.isnan(row).any():
            print(f"WARNING: {drug_name} has NaN in RSV")
            has_warning = True

        if np.all(row == 0.0):
            print(f"WARNING: {drug_name} has all-zero RSV")
            has_warning = True

    return has_warning


def main() -> int:
    repo_root = Path(__file__).parent
    matrix_path = repo_root / "drug_rsv_matrix.npy"
    names_path = repo_root / "drug_names.txt"

    if not matrix_path.exists():
        print(f"ERROR: drug_rsv_matrix.npy not found at {matrix_path}")
        return 1

    matrix = np.load(matrix_path)
    names = _load_names(names_path)

    print("RSV matrix verification")
    print("=======================")
    print(f"Matrix shape: {matrix.shape}")

    nan_rows = np.isnan(matrix).any(axis=1).sum()
    zero_rows = (matrix == 0.0).all(axis=1).sum()

    print(f"Rows with any NaN: {nan_rows}")
    print(f"Rows all zeros: {zero_rows}")

    labels = [
        "Column 0 = Betweenness Shift",
        "Column 1 = Community Change",
        "Column 2 = Spectral Gap Delta",
        "Column 3 = Entropy Delta",
    ]

    _print_col_stats(matrix, labels)
    _print_specific_drugs(matrix, names)

    print("\nWarnings:")
    has_warning = _scan_warnings(matrix, names)

    if not has_warning:
        print("\nMATRIX VERIFIED OK")
    else:
        print("\nMATRIX HAS ISSUES — check warnings above")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
