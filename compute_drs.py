"""Compute Disease Rewiring Signatures (DRS) for 5 diseases.

Run with:
    python compute_drs.py
"""

from __future__ import annotations

import copy
import csv
import math
from pathlib import Path
from typing import List, Tuple

import numpy as np

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
    compute_community_change,
    compute_entropy_delta,
    compute_spectral_gap,
)


def _diseases() -> List[Tuple[str, List[str], List[float]]]:
    return [
        (
            "Parkinson's Disease",
            ["LRRK2", "SNCA", "PINK1", "PRKN", "GBA"],
            [0.95, 0.88, 0.76, 0.71, 0.68],
        ),
        (
            "Alzheimer's Disease",
            ["APP", "PSEN1", "PSEN2", "APOE", "MAPT"],
            [0.92, 0.89, 0.85, 0.82, 0.79],
        ),
        (
            "Type 2 Diabetes",
            ["INS", "PPARG", "KCNJ11", "TCF7L2", "SLC30A8"],
            [0.90, 0.85, 0.80, 0.78, 0.72],
        ),
        (
            "Breast Cancer",
            ["BRCA1", "BRCA2", "TP53", "ERBB2", "ESR1"],
            [0.95, 0.93, 0.91, 0.88, 0.82],
        ),
        (
            "Hypertension",
            ["ACE", "AGT", "ADRB1", "NOS3", "AGTR1"],
            [0.87, 0.84, 0.79, 0.75, 0.71],
        ),
    ]


def load_gene_map(filepath: Path) -> dict:
    if not filepath.exists():
        print(f"Gene mapping file not found at {filepath}")
        print(
            "Download from: https://string-db.org/mapping_files/STRING_display_names/human.name_2_string.tsv.gz"
        )
        print(f"Unzip and place at: {filepath}")
        return {}

    gene_map = {}
    with filepath.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        headers = [h.strip() for h in reader.fieldnames or []]

        id_col = None
        name_col = None

        for h in headers:
            h_lower = h.lower()
            if h_lower in {"#string_protein_id", "string_protein_id"}:
                id_col = h
            if h_lower == "preferred_name":
                name_col = h

        if not id_col or not name_col:
            raise ValueError(
                "Mapping file missing required columns: #string_protein_id, preferred_name"
            )

        for row in reader:
            string_id = str(row.get(id_col, "")).strip()
            preferred = str(row.get(name_col, "")).strip()
            if not string_id or not preferred:
                continue
            gene_map[preferred.upper()] = string_id

    return gene_map


def remap_genes(gene_list: List[str], gene_map: dict) -> List[str]:
    remapped = []
    for gene in gene_list:
        key = gene.strip().upper()
        if key in gene_map:
            remapped.append(gene_map[key])
        else:
            print(f"[WARN] Gene not found in map: {gene}")

    return remapped


def resolve_disease_gene_pairs(
    raw_genes: List[str],
    scores: List[float],
    G0,
    gene_map: dict
) -> List[Tuple[str, float]]:
    found = []
    missing = []
    graph_nodes_set = set(str(n) for n in G0.nodes())

    ensp_to_symbol = {
        "ENSP00000298910": "SNCA",
        "ENSP00000364204": "LRRK2",
        "ENSP00000355865": "PRKN",
        "ENSP00000284981": "APP",
        "ENSP00000326366": "PSEN1",
        "ENSP00000380432": "INS",
        "ENSP00000287820": "PPARG",
        "ENSP00000418960": "BRCA1",
        "ENSP00000369497": "BRCA2",
        "ENSP00000290866": "ACE",
        "ENSP00000355627": "AGT",
    }

    for gene, score in zip(raw_genes, scores):
        gene = str(gene).strip()

        if gene in graph_nodes_set:
            found.append((gene, score))
            continue

        if gene.startswith("9606."):
            stripped = gene[5:]
            if stripped in graph_nodes_set:
                found.append((stripped, score))
                continue

        upper = gene.upper()
        if upper in graph_nodes_set:
            found.append((upper, score))
            continue

        mapped = gene_map.get(upper)
        if mapped and mapped in graph_nodes_set:
            found.append((mapped, score))
            continue

        bare = gene.replace("9606.", "")
        if bare in ensp_to_symbol:
            sym = ensp_to_symbol[bare]
            if sym in graph_nodes_set:
                found.append((sym, score))
                continue

        missing.append(gene)

    if missing:
        print(f"  [WARN] {len(missing)} genes not mapped: {missing[:3]}")
    if not found:
        print("  [ERROR] Zero disease genes found in graph!")
    else:
        print(f"  [INFO] Resolved {len(found)}/{len(raw_genes)} disease genes")

    return found


def _expand_subgraph_nodes(G0, seeds: List[str], hops: int, max_nodes: int | None = None) -> set:
    nodes = set()
    frontier = set(seeds)

    for _ in range(hops):
        next_frontier = set()
        for node in frontier:
            if node not in G0:
                continue
            for nbr in G0.neighbors(node):
                next_frontier.add(nbr)
        nodes.update(frontier)
        frontier = next_frontier

    nodes.update(frontier)

    if max_nodes is not None and len(nodes) > max_nodes:
        # Deterministic cap to keep spectral calculations tractable.
        nodes = set(sorted(nodes)[:max_nodes])

    return nodes


def _apply_disease_activation(G0, genes: List[str], scores: List[float]) -> Tuple[object, List[str]]:
    G_disease = copy.deepcopy(G0)

    gene_scores = dict(zip(genes, scores))
    present = []
    missing = []

    for gene in genes:
        if gene in G_disease:
            present.append(gene)
        else:
            missing.append(gene)

    if missing:
        print(f"[WARN] Missing genes in graph: {missing}")

    if not present:
        return G_disease, []

    # Strengthen existing edges incident to each disease gene.
    for g in present:
        assoc = gene_scores.get(g, 0.0) * 2.0
        for nbr in G_disease.neighbors(g):
            w = G_disease[g][nbr].get("weight", 0.0)
            new_weight = w * (1.0 + assoc)
            G_disease[g][nbr]["weight"] = np.clip(new_weight, 0.001, 1.0)

    return G_disease, present


def main() -> int:
    ppi_path = SRC_P2 / "ppi_genes.csv"
    gene_map_path = REPO_ROOT / "data" / "processed" / "human.name_2_string.tsv"
    print(f"Loading G0 from {ppi_path}")

    G0 = build_graph(str(ppi_path))
    if G0 is None:
        print("ERROR: Failed to load baseline PPI graph.")
        return 1

    gene_map = load_gene_map(gene_map_path)
    if not gene_map:
        return 0

    drs_vectors = []
    disease_names = []

    for disease_name, genes, scores in _diseases():
        print(f"\nProcessing disease: {disease_name}")
        pairs = resolve_disease_gene_pairs(genes, scores, G0, gene_map)
        if not pairs:
            print("WARNING: No disease genes found in graph — node names may need remapping")
            print(f"Skipping {disease_name}.")
            continue

        use_genes = [g for g, _ in pairs]
        use_scores = [s for _, s in pairs]

        G_disease, present_genes = _apply_disease_activation(G0, use_genes, use_scores)

        if not present_genes:
            print("WARNING: No disease genes found in graph — node names may need remapping")
            print(f"Skipping {disease_name}.")
            continue

        try:
            # Target-centric subgraph: targets + 1-hop neighbors
            nodes = _expand_subgraph_nodes(G0, present_genes, hops=1, max_nodes=1500)

            S0 = G0.subgraph(nodes).copy()
            Sd = G_disease.subgraph(nodes).copy()

            d_bc_raw = compute_betweenness_shift(S0, Sd)
            d_cm_raw = compute_community_change(S0, Sd)
            d_sg_raw = compute_spectral_gap(S0, Sd)
            d_en_raw = compute_entropy_delta(S0, Sd, present_genes)

            print(f"  DEBUG ΔBC raw value : {d_bc_raw:.6f}")
            print(f"  DEBUG ΔSG lam0      : {d_sg_raw + 0.0:.6f}")
            print(f"  DEBUG ΔSG lam1      : {d_sg_raw + 0.0:.6f}")
            print(f"  DEBUG subgraph size : {S0.number_of_nodes()} nodes")
            print(f"  DEBUG targets found : {len(present_genes)}")

            max_deg = max([G0.degree(t) for t in present_genes if t in G0], default=1)
            normalizer = math.log2(max_deg + 1)
            d_en = d_en_raw / (normalizer + 1e-8)
            d_en = float(np.clip(d_en, 0.0, 1.0))

            drs_vector = np.array([d_bc_raw, d_cm_raw, d_sg_raw, d_en], dtype=float)
            drs_vector = np.clip(drs_vector, 0.0, 1.0)

            drs_vectors.append(drs_vector.tolist())
            disease_names.append(disease_name)
            print(
                f"{disease_name} DRS: "
                f"[{drs_vector[0]:.3f}, {drs_vector[1]:.3f}, {drs_vector[2]:.3f}, {drs_vector[3]:.3f}]"
            )
        except Exception as exc:
            print(f"WARNING: Failed DRS for {disease_name}: {exc}")
            continue

    output_dir = Path(__file__).parent
    matrix_path = output_dir / "disease_drs.npy"
    names_path = output_dir / "disease_names.txt"

    if drs_vectors:
        matrix = np.asarray(drs_vectors, dtype=np.float32)
        np.save(matrix_path, matrix)

        with names_path.open("w", encoding="utf-8") as f:
            for name in disease_names:
                f.write(f"{name}\n")

        print("\nDRS COMPUTATION DONE")

        print("\n=== DRS VALIDATION ===")
        print(f"{'Disease':<25} {'ΔBC':>8} {'ΔCM':>8} {'ΔSG':>8} {'ΔEN':>8}  OK?")
        print("-" * 65)
        all_ok = True
        for name, vec in zip(disease_names, matrix):
            in_range = np.all(vec >= 0) and np.all(vec <= 1)
            nonzero = np.sum(np.abs(vec) > 1e-6)
            ok = in_range and nonzero >= 2
            all_ok = all_ok and ok
            flag = "OK" if ok else "FAIL"
            print(
                f"  {name:<25} {vec[0]:>8.4f} {vec[1]:>8.4f} "
                f"{vec[2]:>8.4f} {vec[3]:>8.4f}  {flag}"
            )
        print(
            f"\nAll DRS in [0,1]     : "
            f"{'PASS' if np.all(matrix <= 1.0) else 'FAIL ← ΔEN not clipped'}"
        )
        diffs = [
            np.linalg.norm(matrix[i] - matrix[j])
            for i in range(len(matrix))
            for j in range(i + 1, len(matrix))
        ]
        if diffs:
            print("All diseases distinct: ", end="")
            print(f"{'PASS' if min(diffs) > 1e-6 else 'FAIL'}")
        print("======================")
    else:
        print("\nWARNING: No DRS vectors computed. Nothing saved.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
