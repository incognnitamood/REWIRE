#!/usr/bin/env python3
"""
REWIRE Setup Verification Script

Validates that all dependencies are installed and core functionality works.
Run after: pip install -r requirements.txt

Exit codes:
  0 = All checks passed ✓
  1 = Some checks failed ✗
"""

import sys
from pathlib import Path

def check_module(name, import_name=None):
    """Check if a module is installed."""
    if import_name is None:
        import_name = name
    try:
        __import__(import_name)
        print(f"  ✓ {name}")
        return True
    except ImportError as e:
        print(f"  ✗ {name}: {e}")
        return False

def check_file(path, desc):
    """Check if a file exists."""
    if Path(path).exists():
        print(f"  ✓ {desc}")
        return True
    else:
        print(f"  ✗ {desc} (missing: {path})")
        return False

def main():
    print("\n" + "="*60)
    print("REWIRE Setup Verification")
    print("="*60 + "\n")
    
    all_ok = True
    
    # 1. Check core dependencies
    print("1. Checking Python packages...")
    packages = [
        ("numpy", "numpy"),
        ("pandas", "pandas"),
        ("networkx", "networkx"),
        ("scikit-learn", "sklearn"),
        ("scipy", "scipy"),
        ("requests", "requests"),
        ("python-louvain", "community"),
    ]
    
    for name, import_name in packages:
        if not check_module(name, import_name):
            all_ok = False
    
    print()
    
    # 2. Check streamlit (optional)
    print("2. Checking optional packages...")
    if not check_module("streamlit", "streamlit"):
        print("  (Streamlit required for web interface only)")
    print()
    
    # 3. Check core files
    print("3. Checking core files...")
    files = [
        ("rsv_compute.py", "Core RSV metrics"),
        ("setup_test.py", "Test infrastructure"),
        ("requirements.txt", "Dependency list"),
        ("README.md", "User guide"),
        ("SETUP.md", "Setup documentation"),
    ]
    
    root = Path(__file__).parent
    for filename, desc in files:
        if not check_file(root / filename, desc):
            all_ok = False
    
    print()
    
    # 4. Check pre-computed data
    print("4. Checking pre-computed data...")
    data_files = [
        ("drug_rsv_matrix.npy", "Drug RSV matrix"),
        ("disease_drs.npy", "Disease signatures"),
        ("drug_names.txt", "Drug names"),
        ("disease_names.txt", "Disease names"),
    ]
    
    for filename, desc in data_files:
        if not check_file(root / filename, desc):
            all_ok = False
    
    print()
    
    # 5. Check REWIRE submodule
    print("5. Checking REWIRE submodule...")
    rewire_path = root / "REWIRE" / "src" / "P2"
    if check_file(rewire_path / "ppi_graph.py", "REWIRE/src/P2/ppi_graph.py"):
        # Try importing it
        try:
            sys.path.insert(0, str(rewire_path))
            sys.path.insert(0, str(root / "REWIRE" / "src"))
            from ppi_graph import build_graph
            print("  ✓ ppi_graph.build_graph() importable")
        except Exception as e:
            print(f"  ✗ ppi_graph import failed: {e}")
            all_ok = False
    else:
        all_ok = False
    
    print()
    
    # 6. Test core functionality
    print("6. Testing core functionality...")
    try:
        from rsv_compute import compute_rsv
        import networkx as nx
        import random
        random.seed(42)
        
        # Create simple test graphs
        G0 = nx.Graph()
        G0.add_weighted_edges_from([(0, 1, 0.5), (1, 2, 0.8), (2, 3, 0.3)])
        
        G_drug = G0.copy()
        G_drug[0][1]['weight'] = 0.2
        
        result = compute_rsv(G0, G_drug, targets=[0, 1, 2])
        
        if "vector" in result and len(result["vector"]) == 4:
            print(f"  ✓ compute_rsv() works (RSV: {result['vector'][:2]}...)")
        else:
            print("  ✗ compute_rsv() returned invalid structure")
            all_ok = False
            
    except Exception as e:
        print(f"  ✗ compute_rsv() test failed: {e}")
        all_ok = False
    
    print()
    
    # Final report
    if all_ok:
        print("="*60)
        print("✓ ALL CHECKS PASSED - Setup is valid!")
        print("="*60)
        print("\nNext steps:")
        print("  1. python rsv_compute.py           (verify RSV metrics)")
        print("  2. python compute_drs.py           (compute disease signatures)")
        print("  3. streamlit run app.py            (start web interface)")
        print("\nSee SETUP.md for detailed documentation.\n")
        return 0
    else:
        print("="*60)
        print("✗ SOME CHECKS FAILED - See issues above")
        print("="*60)
        print("\nTroubleshooting:")
        print("  1. Run: pip install -r requirements.txt")
        print("  2. Run: pip list | grep -i (module name)")
        print("  3. See SETUP.md for known issues")
        print()
        return 1

if __name__ == "__main__":
    sys.exit(main())
