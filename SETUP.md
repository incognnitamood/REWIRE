# REWIRE Setup Guide

## Quick Start

### 1. Install Dependencies
```bash
pip install -r requirements.txt
```

### 2. Verify Installation
```bash
python rsv_compute.py
```

Expected output:
```
Betweenness shift: 0.0225
Community change: 0.5030
Spectral gap: -0.1556
Entropy delta: -0.0064
Vector: [0.0225, 0.5030, -0.1556, -0.0064]
ALL TESTS PASSED
```

---

## Project Structure

```
REWIRE/
├── requirements.txt              # Python dependencies (INSTALL FIRST)
├── rsv_compute.py                # Core RSV (Repurposing Signature Vector) metrics
├── setup_test.py                 # Test graph generation utility
├── compute_drs.py                # Disease Rewiring Signatures computation
├── disease_drs.py                # Disease analysis (requires OpenTargets API)
├── compute_batch_rsv.py          # Batch RSV computation
├── app.py                         # Streamlit web interface
├── rank_drugs.py                 # Drug ranking/scoring
├── drug_rsv_matrix.npy           # Pre-computed RSV matrix (20 drugs × 4 metrics)
├── disease_drs.npy               # Pre-computed disease signatures
├── REWIRE/                        # Subdirectory with source code
│   ├── src/
│   │   ├── P2/
│   │   │   ├── ppi_graph.py       # PPI graph builder (required by compute_drs.py)
│   │   │   ├── ppi_genes.csv      # Protein interaction data
│   │   │   ├── drug_targets_symbols.csv
│   │   │   └── disease_genes.csv
│   │   └── rewire_rsv.py
│   ├── data/
│   │   └── processed/
│   │       └── canonical_drug_targets.csv
│   └── README.md
└── *.npy, *.txt, *.csv files     # Cached data and results
```

---

## Main Scripts

### A. Core RSV Computation (✅ Works Out-of-Box)
```bash
python rsv_compute.py
```
**Purpose**: Compute 4 topological metrics on protein graphs
- **Betweenness shift**: Centrality redistribution [0, 1]
- **Community change**: Pathway reorganization [0, 1]
- **Spectral gap**: Network resilience (-∞, ∞)
- **Entropy delta**: Neighborhood chaos (-∞, ∞)

**Output**: Dictionary with named metrics + vector form
```python
{
    "betweenness_shift": 0.0225,
    "community_change": 0.5030,
    "spectral_gap": -0.1556,
    "entropy_delta": -0.0064,
    "vector": [0.0225, 0.5030, -0.1556, -0.0064]
}
```

### B. Disease Rewiring Signatures (⚠️ Requires Setup)
```bash
python compute_drs.py
```
**Purpose**: Compute disease signatures for 5 diseases
**Status**: ✅ Ready (uses local PPI data)
**Data**: Uses REWIRE/src/P2/ PPI and disease gene data

### C. Disease Analysis via OpenTargets API (❌ Requires API Key)
```bash
python disease_drs.py
```
**Purpose**: Fetch disease targets from OpenTargets
**Issue**: `400 Bad Request` from OpenTargets API
**Workaround**: Use pre-computed data instead
```bash
# Use pre-computed results
python verify_rsv_matrix.py
python rank_drugs.py
```

### D. Web Interface (⚠️ Requires Streamlit)
```bash
streamlit run app.py
```
**Purpose**: Interactive drug-disease matching
**Data**: Uses pre-computed matrices (drug_rsv_matrix.npy, disease_drs.npy)

---

## Data Files

Pre-computed files (included):
- `drug_rsv_matrix.npy` - RSV vectors for 20 drugs
- `disease_drs.npy` - Disease signatures for 5 diseases
- `drug_names.txt` - Drug identifiers
- `disease_names.txt` - Disease names
- `chembl_targets_raw.csv` - ChEMBL target data

These allow you to run analysis without recomputing.

---

## Known Issues & Solutions

### ❌ Issue 1: `ModuleNotFoundError: No module named 'ppi_graph'`
**Cause**: compute_drs.py imports from REWIRE/src/P2/ but Python path isn't set
**Solution** (ALREADY IN CODE): 
```python
sys.path.insert(0, str(SRC_P2))
```
✅ compute_drs.py already has this. If you get this error:
1. Verify REWIRE/src/P2/ppi_graph.py exists
2. Run from root directory (c:\Users\sujat\OneDrive\Desktop\REWIRE\)

### ❌ Issue 2: OpenTargets API 400 Bad Request
**Cause**: API endpoint changed or authentication needed
**Workaround**: Use pre-computed disease_drs.npy instead
```python
# Instead of disease_drs.py, use:
import numpy as np
disease_drs = np.load("disease_drs.npy")
```

### ❌ Issue 3: Missing `ppi_graph` module in rsv_compute.py
**Why**: rsv_compute.py is standalone (no external imports except networkx, scipy, sklearn)
✅ **Status**: Fixed - rsv_compute.py is production-ready

### ⚠️ Issue 4: Requirements.txt was empty
**Why**: Hadn't been populated
✅ **Status**: Fixed - now contains all dependencies

---

## Installation Verification

After `pip install -r requirements.txt`, verify all modules:

```bash
python -c "import networkx, numpy, pandas, sklearn, scipy, streamlit, community; print('✓ All core modules installed')"
```

---

## Running Tests

```bash
# Test RSV computation (core metrics)
python rsv_compute.py

# Test graph setup
python setup_test.py

# Test RSV matrix verification
python verify_rsv_matrix.py

# Test single RSV computation
python test_single_pkl.py
```

---

## For Developers

### Adding New Dependencies
1. Update `requirements.txt`
2. Run `pip install -r requirements.txt`
3. Test with: `python rsv_compute.py`

### Running Disease Analysis Locally
```python
from REWIRE.src.P2.ppi_graph import build_graph
from rsv_compute import compute_rsv

# Build PPI graph
G = build_graph()

# Compute RSV for a disease
rsv = compute_rsv(G, G_modified, targets=["LRRK2", "SNCA"])
print(rsv["vector"])
```

### Using Pre-computed Data
```python
import numpy as np

# Load drug RSV matrix (20 drugs × 4 metrics)
drug_rsv = np.load("drug_rsv_matrix.npy")

# Load disease signatures (5 diseases)
disease_drs = np.load("disease_drs.npy")

# Compute similarity scores
from sklearn.metrics.pairwise import cosine_similarity
scores = cosine_similarity(drug_rsv, disease_drs)
```

---

## API Documentation

### `compute_rsv(G0, G_drug, targets) -> dict`

Returns:
```python
{
    "betweenness_shift": float,   # [0, 1] - hub centrality change
    "community_change": float,     # [0, 1] - module reorganization
    "spectral_gap": float,         # unbounded - network fragmentation
    "entropy_delta": float,        # unbounded - target chaos
    "vector": [4 floats]           # For ML pipelines
}
```

**On Failure**: Component = `NaN` (detectable, not silent 0.0)

---

## Contact & Issues

- **rsv_compute.py**: Fully tested ✅
- **compute_drs.py**: Tested with local data ✅
- **disease_drs.py**: Requires external API setup ⚠️
- **app.py**: Web interface (use pre-computed data) ⚠️

For module/import issues, verify REWIRE/ subdirectory structure matches the paths in sys.path.insert() calls.
