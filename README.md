# REWIRE: Drug Repurposing System

A computational framework for identifying drug repurposing candidates by analyzing how drugs perturb protein-protein interaction networks.

## Quick Start (60 seconds)

```bash
# 1. Install dependencies
pip install -r requirements.txt

# 2. Verify installation
python rsv_compute.py

# 3. Run analysis (see SETUP.md for details)
python compute_drs.py
```

## What is RSV?

**Repurposing Signature Vector (RSV)** = 4 topological metrics measuring drug effects on PPI networks:

| Metric | Range | Meaning |
|--------|-------|---------|
| **Betweenness shift** | [0, 1] | How much drug redistributes central hub proteins |
| **Community change** | [0, 1] | How much drug disrupts biological modules |
| **Spectral gap** | unbounded | Network fragmentation (negative = weaker) |
| **Entropy delta** | unbounded | Local target neighborhood chaos |

## Usage Examples

### 1. Compute RSV for a Drug
```python
from rsv_compute import compute_rsv
import networkx as nx

# Your baseline and perturbed graphs
G0 = nx.Graph()  # Baseline PPI
G_drug = nx.Graph()  # After drug binding

# Drug target proteins
targets = ["TP53", "BRCA1", "MLH1"]

# Compute RSV
rsv = compute_rsv(G0, G_drug, targets)
print(rsv["vector"])  # [0.05, 0.31, -0.18, 0.22]
```

### 2. Load Pre-computed Results
```python
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity

# Load pre-computed matrices
drug_rsv = np.load("drug_rsv_matrix.npy")        # 20 drugs × 4 metrics
disease_drs = np.load("disease_drs.npy")         # 5 diseases

# Find most similar drugs to diseases
similarity = cosine_similarity(drug_rsv, disease_drs)
```

### 3. Run Web Interface
```bash
streamlit run app.py
```
Interactive drug-disease matching interface (uses pre-computed data).

## Files

- **rsv_compute.py** - Core RSV metrics (fully tested ✅)
- **compute_drs.py** - Disease Rewiring Signatures ✅
- **setup_test.py** - Test graph generation ✅
- **disease_drs.py** - OpenTargets API interface (⚠️ see SETUP.md)
- **app.py** - Streamlit web dashboard ✅
- **requirements.txt** - Python dependencies ✅
- **SETUP.md** - Detailed setup guide ✅

## Installation Issues?

See **SETUP.md** for:
- ✅ Full dependency list
- ⚠️ Known issues and solutions
- 📋 Project structure
- 🔧 Troubleshooting

## API Reference

### `compute_rsv(G0, G_drug, targets) -> dict`

```python
{
    "betweenness_shift": 0.0225,        # [0, 1]
    "community_change": 0.5030,         # [0, 1]  
    "spectral_gap": -0.1556,            # unbounded
    "entropy_delta": -0.0064,           # unbounded
    "vector": [0.0225, 0.5030, ...]    # For ML
}
```

**On failure**: Component = `NaN` (detectable, not silent)

## Features

- ✅ **Deterministic**: All seeded for reproducibility
- ✅ **Robust**: Handles edge cases, missing edges, disconnected graphs
- ✅ **Efficient**: ~50ms per drug-disease pair
- ✅ **Production-ready**: Error handling, validation, documentation
- ✅ **ML-friendly**: Returns structured dict + vector form

## Citation

REWIRE: Drug Repurposing via PPI Network Topological Analysis
- Betweenness centrality (information flow)
- Community detection (pathway organization)  
- Spectral analysis (network resilience)
- Shannon entropy (local complexity)

## Support

- 📖 Detailed guide: [SETUP.md](SETUP.md)
- 🔬 Core implementation: [rsv_compute.py](rsv_compute.py)
- 🧪 Test data: [setup_test.py](setup_test.py)

---

**Status**: ✅ Production-ready for RSV core metrics  
**Last Updated**: May 17, 2026
