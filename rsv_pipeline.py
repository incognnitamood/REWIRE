#!/usr/bin/env python3
"""
REWIRE RSV Pipeline: Unified workflow for drug repurposing analysis

This script orchestrates the complete analysis pipeline:
1. Compute Disease Rewiring Signatures (DRS) from OpenTargets API
2. Compute Drug RSV vectors from PPI network perturbations
3. Rank drugs by similarity to disease signatures

Usage:
    python rsv_pipeline.py [--diseases DISEASE1,DISEASE2,... | --pre-computed]
    python rsv_pipeline.py --help
    
Examples:
    # Use pre-computed disease signatures
    python rsv_pipeline.py --pre-computed
    
    # Fetch fresh disease data from OpenTargets API
    python rsv_pipeline.py --diseases "Parkinson's Disease,Alzheimer's Disease"
    
    # Compute all steps from scratch
    python rsv_pipeline.py --compute-all
"""

import argparse
import sys
from pathlib import Path
from typing import Optional

import numpy as np

# Ensure imports work
REPO_ROOT = Path(__file__).parent / "REWIRE"
SRC_ROOT = REPO_ROOT / "src"
SRC_P2 = SRC_ROOT / "P2"
sys.path.insert(0, str(SRC_ROOT))
sys.path.insert(0, str(SRC_P2))


def load_or_compute_disease_drs(
    diseases: Optional[list] = None,
    force_recompute: bool = False
) -> tuple:
    """
    Load disease signatures from cache or compute fresh from OpenTargets.
    
    Args:
        diseases: List of disease names to compute. If None, uses pre-computed.
        force_recompute: Force recomputation even if cache exists.
    
    Returns:
        (disease_drs: np.ndarray shape (n_diseases, 4), disease_names: list)
    """
    cache_path = Path("disease_drs.npy")
    names_path = Path("disease_names.txt")
    
    # If using pre-computed and cache exists
    if not force_recompute and diseases is None:
        if cache_path.exists() and names_path.exists():
            print("✓ Loading pre-computed disease signatures from cache...")
            disease_drs = np.load(str(cache_path))
            with open(str(names_path)) as f:
                disease_names = [line.strip() for line in f.readlines() if line.strip()]
            print(f"  Loaded {len(disease_names)} disease signatures")
            return disease_drs, disease_names
        else:
            print("⚠ Pre-computed disease signatures not found!")
            print("  Use: python rsv_pipeline.py --diseases DISEASE1,DISEASE2,...")
            sys.exit(1)
    
    # Compute fresh disease signatures
    if diseases is None:
        diseases = [
            "Parkinson's Disease",
            "Alzheimer's Disease",
            "Type 2 Diabetes",
            "Breast Cancer",
            "Hypertension",
        ]
    
    print(f"\n📊 Computing Disease Rewiring Signatures ({len(diseases)} diseases)...")
    
    try:
        from disease_drs import compute_drs
    except Exception as e:
        print(f"✗ Failed to import disease_drs: {e}")
        print("  Falling back to pre-computed data...")
        if cache_path.exists() and names_path.exists():
            disease_drs = np.load(str(cache_path))
            with open(str(names_path)) as f:
                disease_names = [line.strip() for line in f.readlines() if line.strip()]
            return disease_drs, disease_names
        sys.exit(1)
    
    drs_vectors = []
    valid_diseases = []
    
    for disease in diseases:
        try:
            print(f"  Processing: {disease}...", end=" ", flush=True)
            drs_vector = compute_drs(disease)
            if np.allclose(drs_vector, 0.0):
                print("⚠ (zero vector - API may have failed)")
            else:
                print(f"✓ ({drs_vector[0]:.4f}, {drs_vector[1]:.4f}, ...)")
            drs_vectors.append(drs_vector)
            valid_diseases.append(disease)
        except Exception as e:
            print(f"✗ Error: {e}")
    
    if not drs_vectors:
        print("✗ No disease signatures computed!")
        sys.exit(1)
    
    disease_drs = np.vstack(drs_vectors)
    
    # Save for future use
    np.save(str(cache_path), disease_drs)
    with open(str(names_path), "w") as f:
        f.write("\n".join(valid_diseases))
    
    print(f"✓ Computed {len(valid_diseases)} disease signatures")
    print(f"  Saved to: {cache_path}, {names_path}")
    
    return disease_drs, valid_diseases


def load_or_compute_drug_rsv(force_recompute: bool = False) -> tuple:
    """
    Load drug RSV matrix from cache or compute fresh.
    
    Args:
        force_recompute: Force recomputation even if cache exists.
    
    Returns:
        (drug_rsv: np.ndarray shape (n_drugs, 4), drug_names: list)
    """
    cache_path = Path("drug_rsv_matrix.npy")
    names_path = Path("drug_names.txt")
    
    # Check for pre-computed cache
    if not force_recompute and cache_path.exists() and names_path.exists():
        print("\n✓ Loading pre-computed drug RSV matrix...")
        drug_rsv = np.load(str(cache_path))
        with open(str(names_path)) as f:
            drug_names = [line.strip() for line in f.readlines() if line.strip()]
        print(f"  Loaded {len(drug_names)} drug RSV vectors ({drug_rsv.shape})")
        return drug_rsv, drug_names
    
    # Try to compute fresh
    print("\n🔬 Computing Drug RSV Matrix...")
    try:
        from compute_batch_rsv import compute_batch_rsv
        drug_rsv, drug_names = compute_batch_rsv()
        np.save(str(cache_path), drug_rsv)
        with open(str(names_path), "w") as f:
            f.write("\n".join(drug_names))
        print(f"✓ Computed {len(drug_names)} drug RSV vectors")
        return drug_rsv, drug_names
    except Exception as e:
        print(f"⚠ Could not compute drug RSVs: {e}")
        print("  Falling back to pre-computed matrix...")
        if cache_path.exists() and names_path.exists():
            drug_rsv = np.load(str(cache_path))
            with open(str(names_path)) as f:
                drug_names = [line.strip() for line in f.readlines() if line.strip()]
            return drug_rsv, drug_names
        print("✗ No drug RSV matrix available!")
        sys.exit(1)


def rank_drugs(
    drug_rsv: np.ndarray,
    disease_drs: np.ndarray,
    drug_names: list,
    disease_names: list,
    top_k: int = 10
) -> dict:
    """
    Rank drugs by similarity to disease signatures.
    
    Args:
        drug_rsv: Drug RSV matrix (n_drugs, 4)
        disease_drs: Disease DRS matrix (n_diseases, 4)
        drug_names: List of drug names
        disease_names: List of disease names
        top_k: Number of top drugs to return per disease
    
    Returns:
        Dictionary mapping disease names to ranked drug lists
    """
    from sklearn.metrics.pairwise import cosine_similarity
    
    print("\n🎯 Ranking drugs by disease similarity...")
    
    # Compute cosine similarity
    similarity = cosine_similarity(drug_rsv, disease_drs)  # (n_drugs, n_diseases)
    
    results = {}
    for disease_idx, disease_name in enumerate(disease_names):
        scores = similarity[:, disease_idx]
        top_indices = np.argsort(-scores)[:top_k]  # Sort descending
        
        top_drugs = [
            {
                "rank": i + 1,
                "drug": drug_names[idx],
                "score": float(scores[idx]),
            }
            for i, idx in enumerate(top_indices)
        ]
        
        results[disease_name] = top_drugs
        
        print(f"\n  {disease_name}:")
        for drug_info in top_drugs:
            print(f"    {drug_info['rank']}. {drug_info['drug']}: {drug_info['score']:.4f}")
    
    return results


def main():
    parser = argparse.ArgumentParser(
        description="REWIRE RSV Pipeline - Drug repurposing analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python rsv_pipeline.py --pre-computed
      Use pre-computed disease signatures and drug RSVs
  
  python rsv_pipeline.py --diseases "Parkinson's Disease,Alzheimer's Disease"
      Fetch fresh disease data from OpenTargets API
  
  python rsv_pipeline.py --compute-all
      Recompute everything from scratch
        """,
    )
    
    parser.add_argument(
        "--pre-computed",
        action="store_true",
        help="Use pre-computed disease signatures and drug RSVs (default)"
    )
    parser.add_argument(
        "--diseases",
        type=str,
        help='Disease names (comma-separated) to fetch from OpenTargets API'
    )
    parser.add_argument(
        "--compute-all",
        action="store_true",
        help="Force recomputation of all signatures"
    )
    parser.add_argument(
        "--top-k",
        type=int,
        default=10,
        help="Number of top drugs to return per disease (default: 10)"
    )
    
    args = parser.parse_args()
    
    print("=" * 70)
    print("REWIRE: Drug Repurposing via PPI Network Analysis")
    print("=" * 70)
    
    # Step 1: Load or compute disease signatures
    if args.diseases:
        diseases = [d.strip() for d in args.diseases.split(",")]
        disease_drs, disease_names = load_or_compute_disease_drs(
            diseases=diseases,
            force_recompute=args.compute_all
        )
    else:
        disease_drs, disease_names = load_or_compute_disease_drs(
            force_recompute=args.compute_all
        )
    
    # Step 2: Load or compute drug RSV matrix
    drug_rsv, drug_names = load_or_compute_drug_rsv(
        force_recompute=args.compute_all
    )
    
    # Step 3: Rank drugs by similarity
    rankings = rank_drugs(
        drug_rsv,
        disease_drs,
        drug_names,
        disease_names,
        top_k=args.top_k
    )
    
    # Save results
    print("\n" + "=" * 70)
    print("✓ Pipeline complete!")
    print("=" * 70)
    print(f"\nSummary:")
    print(f"  Diseases: {len(disease_names)}")
    print(f"  Drugs: {len(drug_names)}")
    print(f"  Similarity metric: Cosine similarity (drug RSV vs disease DRS)")
    
    # Optionally save rankings to file
    try:
        import json
        with open("rankings.json", "w") as f:
            json.dump(rankings, f, indent=2)
        print(f"\nRankings saved to: rankings.json")
    except Exception as e:
        print(f"\n⚠ Could not save rankings.json: {e}")
    
    print("\nNext steps:")
    print("  1. Review rankings above")
    print("  2. Validate top candidates with literature/databases")
    print("  3. Proceed to wet lab validation for promising candidates")


if __name__ == "__main__":
    main()
