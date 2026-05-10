"""Diagnose RSV computation issues by inspecting graphs and metrics.

Run with:
    python diagnose_rsv.py
"""

from __future__ import annotations

import pickle
from pathlib import Path
from typing import Dict, Iterable, Tuple

import numpy as np
import networkx as nx
import pandas as pd

# Add repo src paths for imports
import sys

REPO_ROOT = Path(__file__).parent / "REWIRE"
SRC_ROOT = REPO_ROOT / "src"
SRC_P2 = SRC_ROOT / "P2"

sys.path.insert(0, str(SRC_ROOT))
sys.path.insert(0, str(SRC_P2))

from ppi_graph import build_graph
from rsv_compute import (
    compute_betweenness_shift,
    compute_spectral_gap,
    compute_entropy_delta,
)


def _first_pkl(path: Path) -> Path | None:
    files = sorted(path.glob("*.pkl"))
    return files[0] if files else None


def _print_edge_samples(G: nx.Graph, label: str) -> None:
    print(f"\n{label} edge samples:")
    samples = list(G.edges(data=True))[:3]
    for u, v, data in samples:
        print(u, v, data)


def _find_numeric_attr(G: nx.Graph) -> str | None:
    for _, _, data in G.edges(data=True):
        for k, v in data.items():
            if isinstance(v, (int, float)):
                if k == "weight":
                    return "weight"
        for k, v in data.items():
            if isinstance(v, (int, float)):
                return k
    return None


def _min_max_attr(G: nx.Graph, attr: str) -> Tuple[float, float]:
    values = [data.get(attr) for _, _, data in G.edges(data=True) if isinstance(data.get(attr), (int, float))]
    if not values:
        return (float("nan"), float("nan"))
    return (min(values), max(values))


def _print_neighbors(G: nx.Graph, node: str, attr: str) -> None:
    neighbors = list(G.neighbors(node))
    print(f"Neighbors for {node} ({len(neighbors)}):")
    for nbr in neighbors:
        data = G.get_edge_data(node, nbr) or {}
        print(f"  {node} -- {nbr}  {attr}={data.get(attr)}")


def _entropy_for_target(G: nx.Graph, target: str, attr: str) -> float:
    neighbors = list(G.neighbors(target))
    if not neighbors:
        return 0.0

    weights = []
    for nbr in neighbors:
        data = G.get_edge_data(target, nbr) or {}
        w = data.get(attr, 0.0)
        if not isinstance(w, (int, float)):
            w = 0.0
        weights.append(w)

    total = sum(weights)
    if total == 0:
        return 0.0

    probs = [w / total for w in weights]
    return float(-sum(p * np.log(p) for p in probs if p > 0))


def main() -> int:
    ppi_path = SRC_P2 / "ppi_genes.csv"
    drug_graphs_dir = SRC_P2 / "drug_graphs"
    drug_50_path = SRC_P2 / "drug_50.csv"

    print(f"Loading G0 from {ppi_path}")
    G0 = build_graph(str(ppi_path))
    if G0 is None:
        print("ERROR: Failed to load G0")
        return 1

    first_pkl = _first_pkl(drug_graphs_dir)
    if first_pkl is None:
        print(f"ERROR: No .pkl files found in {drug_graphs_dir}")
        return 1

    print(f"Loading first PKL: {first_pkl}")
    with first_pkl.open("rb") as f:
        G_drug = pickle.load(f)

    _print_edge_samples(G0, "G0")
    _print_edge_samples(G_drug, "G_drug")

    g0_attr = _find_numeric_attr(G0)
    gd_attr = _find_numeric_attr(G_drug)

    print(f"\nG0 numeric edge attribute: {g0_attr}")
    print(f"G_drug numeric edge attribute: {gd_attr}")

    if g0_attr:
        g0_min, g0_max = _min_max_attr(G0, g0_attr)
        print(f"G0 {g0_attr} min={g0_min} max={g0_max}")

    if gd_attr:
        gd_min, gd_max = _min_max_attr(G_drug, gd_attr)
        print(f"G_drug {gd_attr} min={gd_min} max={gd_max}")

    # Node set comparison
    g0_nodes = set(G0.nodes())
    gd_nodes = set(G_drug.nodes())
    print(f"\nNodes in G0 but not G_drug: {len(g0_nodes - gd_nodes)}")
    print(f"Nodes in G_drug but not G0: {len(gd_nodes - g0_nodes)}")

    # Target protein from drug_50.csv
    if not drug_50_path.exists():
        print(f"ERROR: Missing drug_50.csv at {drug_50_path}")
        return 1

    df = pd.read_csv(drug_50_path)
    if df.empty or "target_gene" not in df.columns:
        print("ERROR: drug_50.csv missing target_gene column")
        return 1

    target = str(df.iloc[0]["target_gene"]).strip()
    print(f"\nTarget protein from drug_50.csv: {target}")
    print(f"Target in G0: {target in G0}")
    print(f"Target in G_drug: {target in G_drug}")

    if g0_attr and target in G0:
        print("\nG0 neighbor weights:")
        _print_neighbors(G0, target, g0_attr)
    if gd_attr and target in G_drug:
        print("\nG_drug neighbor weights:")
        _print_neighbors(G_drug, target, gd_attr)

    if g0_attr and gd_attr and target in G0 and target in G_drug:
        changed = False
        neighbors = set(G0.neighbors(target)).union(set(G_drug.neighbors(target)))
        for nbr in neighbors:
            w0 = (G0.get_edge_data(target, nbr) or {}).get(g0_attr, 0.0)
            wd = (G_drug.get_edge_data(target, nbr) or {}).get(gd_attr, 0.0)
            if isinstance(w0, (int, float)) and isinstance(wd, (int, float)) and w0 != wd:
                changed = True
                break
        print(f"\nAny weights changed for target {target}: {changed}")

    # Manual RSV metrics
    print("\nRSV metric diagnostics:")
    try:
        bet = compute_betweenness_shift(G0, G_drug)
        print(f"compute_betweenness_shift: {bet}")
    except Exception as exc:
        print(f"compute_betweenness_shift ERROR: {exc}")

    try:
        ac_g0 = nx.algebraic_connectivity(G0, weight=g0_attr or "weight")
        ac_gd = nx.algebraic_connectivity(G_drug, weight=gd_attr or "weight")
        print(f"algebraic_connectivity(G0): {ac_g0}")
        print(f"algebraic_connectivity(G_drug): {ac_gd}")
        print(f"compute_spectral_gap: {compute_spectral_gap(G0, G_drug)}")
    except Exception as exc:
        print(f"compute_spectral_gap ERROR: {exc}")

    try:
        entropy_delta = compute_entropy_delta(G0, G_drug, [target])
        print(f"compute_entropy_delta (target={target}): {entropy_delta}")
        if g0_attr and target in G0:
            print(f"entropy(G0) for {target}: {_entropy_for_target(G0, target, g0_attr)}")
        if gd_attr and target in G_drug:
            print(f"entropy(G_drug) for {target}: {_entropy_for_target(G_drug, target, gd_attr)}")
    except Exception as exc:
        print(f"compute_entropy_delta ERROR: {exc}")

    print("\nDIAGNOSIS COMPLETE")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
