import numpy as np
import pickle
import pandas as pd
from pathlib import Path
from tqdm import tqdm
import sys, os

# ── adjust these paths if needed ──────────────────────────────────────────────
G0_PATH      = r"C:\Users\sujat\OneDrive\Desktop\REWIRE\REWIRE\src\P2\ppi_genes.csv"
DRUG_CSV     = r"C:\Users\sujat\OneDrive\Desktop\REWIRE\REWIRE\src\P2\drug_50.csv"
MATRIX_PATH  = r"C:\Users\sujat\OneDrive\Desktop\REWIRE\drug_rsv_matrix.npy"
IDS_PATH     = r"C:\Users\sujat\OneDrive\Desktop\REWIRE\drug_ids.txt"
NAMES_PATH   = r"C:\Users\sujat\OneDrive\Desktop\REWIRE\drug_names.txt"
# ──────────────────────────────────────────────────────────────────────────────

sys.path.insert(0, r"C:\Users\sujat\OneDrive\Desktop\REWIRE\REWIRE\src\P2")
from ppi_graph import build_graph, simulate_binding
from rsv_compute import compute_spectral_gap, compute_entropy_delta

print("Loading G0...")
G0 = build_graph(G0_PATH)

matrix = np.load(MATRIX_PATH)
drug_ids   = open(IDS_PATH).read().strip().splitlines()
drug_names = open(NAMES_PATH).read().strip().splitlines()
df = pd.read_csv(DRUG_CSV)

# build lookup: drug_name -> list of targets
# CSV has one row per target (drug_name, target_gene)
lookup = {}
for _, row in df.iterrows():
    name = str(row['drug_name']).strip()
    gene = str(row['target_gene']).strip()
    if name not in lookup:
        lookup[name] = []
    lookup[name].append(gene)

# no affinity values in CSV — use default 1.0 nM for all targets
# inhibition_factor = 1/(1+1) = 0.5, so weights halved (matches what we saw)

print(f"Patching columns 2 and 3 for {len(drug_names)} drugs...")
failed = []

for i, name in enumerate(tqdm(drug_names)):
    if name not in lookup:
        print(f"WARNING: {name} not in drug CSV — skipping")
        failed.append(name)
        continue
    targets = lookup[name]
    affinities = [1.0] * len(targets)
    try:
        affinity_dict = {t: a for t, a in zip(targets, affinities)}
        result = simulate_binding(G0, targets, affinity_dict)
        G_drug = result[0] if isinstance(result, tuple) else result
        t = targets[0]
        if t in G0 and t in G_drug:
            neighbors = list(G0.neighbors(t))
            if neighbors:
                n = neighbors[0]
                w0 = G0[t][n]['weight']
                wd = G_drug[t][n]['weight']
                if w0 == wd:
                    print(f"WARNING {name}: weights unchanged after simulation")
                else:
                    print(f"OK {name}: {w0:.3f} -> {wd:.3f}")
        sg = compute_spectral_gap(G0, G_drug)
        en = compute_entropy_delta(G0, G_drug, targets)
        matrix[i, 2] = sg
        matrix[i, 3] = en
    except Exception as e:
        print(f"WARNING: {name} failed — {e}")
        failed.append(name)

np.save(MATRIX_PATH, matrix)
print(f"\nDone. Failed: {failed}")
print(f"Matrix shape: {matrix.shape}")

# sanity print
for check in ['Imatinib', 'Nilotinib', 'Metformin', 'Bexarotene']:
    if check in drug_names:
        idx = drug_names.index(check)
        print(f"{check} RSV: {matrix[idx].tolist()}")
