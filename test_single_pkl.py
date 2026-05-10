"""
Test script to validate single PKL loading and RSV computation.

This script:
1. Loads the baseline PPI graph (G0) from ppi_graph.py
2. Reads the first drug ID from drug_ids_completed.txt
3. Loads the corresponding pickled drug-perturbed graph
4. Extracts target proteins from the embedded metadata
5. Computes RSV using rsv_compute.py
6. Prints the drug name and RSV vector

Usage:
    python test_single_pkl.py
"""

import sys
import pickle
from pathlib import Path

# Add src to path so we can import modules
sys.path.insert(0, str(Path(__file__).parent / "REWIRE" / "src"))
sys.path.insert(0, str(Path(__file__).parent / "REWIRE" / "src" / "P2"))

# Import required modules
from ppi_graph import build_graph
from rsv_compute import compute_rsv
import pandas as pd


def main():
    """Main test execution."""
    
    print("\n" + "="*70)
    print("SINGLE PKL TEST — Load drug graph and compute RSV")
    print("="*70 + "\n")
    
    # ─────────────────────────────────────────────────────────────────────
    # 1. Load baseline PPI graph (G0)
    # ─────────────────────────────────────────────────────────────────────
    
    repo_root = Path(__file__).parent / "REWIRE"
    ppi_path = repo_root / "src" / "P2" / "ppi_genes.csv"
    
    print(f"[1/6] Loading PPI graph from: {ppi_path}")
    if not ppi_path.exists():
        print(f"ERROR: PPI file not found at {ppi_path}")
        return False
    
    G0 = build_graph(str(ppi_path))
    if G0 is None:
        print("ERROR: Failed to load G0")
        return False
    
    print(f"      ✓ G0 loaded: {G0.number_of_nodes():,} nodes, {G0.number_of_edges():,} edges\n")
    
    # ─────────────────────────────────────────────────────────────────────
    # 2. Read first drug ID from drug_ids_completed.txt
    # ─────────────────────────────────────────────────────────────────────
    
    drug_ids_file = repo_root / "src" / "P2" / "drug_graphs" / "drug_ids_completed.txt"
    
    print(f"[2/6] Reading first drug from: {drug_ids_file}")
    if not drug_ids_file.exists():
        print(f"ERROR: Drug IDs file not found at {drug_ids_file}")
        return False
    
    with open(drug_ids_file, 'r') as f:
        first_drug_name = f.readline().strip()
    
    if not first_drug_name:
        print("ERROR: No drugs found in drug_ids_completed.txt")
        return False
    
    print(f"      ✓ First drug: {first_drug_name}\n")
    
    # ─────────────────────────────────────────────────────────────────────
    # 3. Load corresponding PKL file
    # ─────────────────────────────────────────────────────────────────────
    
    pkl_filename = first_drug_name.lower().replace(" ", "_") + ".pkl"
    pkl_path = repo_root / "src" / "P2" / "drug_graphs" / pkl_filename
    
    print(f"[3/6] Loading drug graph from: {pkl_path}")
    if not pkl_path.exists():
        # Try alternative naming (no replacement)
        pkl_path = repo_root / "src" / "P2" / "drug_graphs" / (first_drug_name.lower() + ".pkl")
        if not pkl_path.exists():
            print(f"ERROR: PKL file not found at {pkl_path}")
            return False
    
    try:
        with open(pkl_path, 'rb') as f:
            G_drug = pickle.load(f)
        print(f"      ✓ G_drug loaded: {G_drug.number_of_nodes():,} nodes, {G_drug.number_of_edges():,} edges\n")
    except Exception as e:
        print(f"ERROR: Failed to load PKL file: {e}")
        return False
    
    # ─────────────────────────────────────────────────────────────────────
    # 4. Extract target proteins from embedded metadata
    # ─────────────────────────────────────────────────────────────────────
    
    print(f"[4/6] Extracting target proteins from metadata")
    
    targets = []
    if hasattr(G_drug, 'graph') and 'rewire' in G_drug.graph:
        metadata = G_drug.graph['rewire']
        if 'targets' in metadata:
            targets = metadata['targets']
            if isinstance(targets, str):
                targets = [t.strip() for t in targets.split(',')]
    
    if not targets:
        print("      WARNING: No targets found in metadata or metadata missing")
        print("      Using empty target list (will compute RSV on full graphs)")
        targets = []
    else:
        print(f"      ✓ Targets found: {targets}\n")
    
    # ─────────────────────────────────────────────────────────────────────
    # 5. Validate targets are in G0 (warn if missing)
    # ─────────────────────────────────────────────────────────────────────
    
    print(f"[5/6] Validating targets against G0")
    
    missing_targets = []
    for target in targets:
        if target not in G0.nodes():
            missing_targets.append(target)
    
    if missing_targets:
        print(f"      WARNING: {len(missing_targets)} targets not found in G0:")
        for t in missing_targets[:5]:  # Print first 5
            print(f"        - {t}")
        if len(missing_targets) > 5:
            print(f"        ... and {len(missing_targets) - 5} more")
        valid_targets = [t for t in targets if t in G0.nodes()]
        if valid_targets:
            print(f"      Using {len(valid_targets)} valid targets: {valid_targets}")
            targets = valid_targets
        else:
            print(f"      No valid targets found. Computing RSV on full graphs.")
            targets = []
    else:
        print(f"      ✓ All {len(targets)} targets found in G0\n")
    
    # ─────────────────────────────────────────────────────────────────────
    # 6. Compute RSV
    # ─────────────────────────────────────────────────────────────────────
    
    print(f"[6/6] Computing RSV for {first_drug_name}")
    print(f"      (This may take 30-120 seconds for betweenness/community detection...)\n")
    
    try:
        rsv_dict = compute_rsv(G0, G_drug, targets)
        vector = rsv_dict.get("vector", [])
        
        print("      ✓ RSV computation succeeded!\n")
        print(f"{'='*70}")
        print(f"RESULT: {first_drug_name}")
        print(f"{'='*70}")
        print(f"  Betweenness Shift:     {rsv_dict['betweenness_shift']:.6f}")
        print(f"  Community Change:      {rsv_dict['community_change']:.6f}")
        print(f"  Spectral Gap:          {rsv_dict['spectral_gap']:.6f}")
        print(f"  Entropy Delta:         {rsv_dict['entropy_delta']:.6f}")
        print(f"\n  RSV Vector: {vector}")
        print(f"{'='*70}\n")
        
        # Validate vector has 4 elements and no NaN
        if len(vector) != 4:
            print(f"ERROR: Vector has {len(vector)} elements, expected 4")
            return False
        
        has_nan = any(v != v for v in vector)  # NaN != NaN
        if has_nan:
            print(f"WARNING: Vector contains NaN values")
            print(f"This may indicate issues with graph structure or connectivity.")
        
        print("✅ SINGLE PKL TEST PASSED\n")
        return True
        
    except Exception as e:
        print(f"ERROR: RSV computation failed: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
