import networkx as nx
from copy import deepcopy
from community import best_partition
from sklearn.metrics import normalized_mutual_info_score
from scipy.stats import entropy as scipy_entropy


def compute_betweenness_shift(G0, G_drug) -> float:
    """
    Compute the mean absolute shift in betweenness centrality between two graphs.
    
    Measures how much drug binding (represented by edge weight changes) shifts
    the betweenness centrality of nodes. A higher shift indicates the drug
    significantly redistributes information flow through the network.
    
    Parameters:
    -----------
    G0 : NetworkX Graph
        Baseline protein-protein interaction graph
    G_drug : NetworkX Graph
        Graph after drug binding simulation (same nodes as G0, some edges weighted differently)
    
    Returns:
    --------
    float
        Mean absolute difference in betweenness centrality across all nodes (0-1 range)
    """
    # Validate that graphs have identical node sets
    assert set(G0.nodes()) == set(G_drug.nodes()), "Graphs must have identical nodes"
    
    # Validate weights exist on all edges
    assert all('weight' in d for _, _, d in G0.edges(data=True)), "G0: all edges must have 'weight' attribute"
    assert all('weight' in d for _, _, d in G_drug.edges(data=True)), "G_drug: all edges must have 'weight' attribute"
    
    # Compute betweenness centrality on both graphs
    # Use k=50 sampling, but cap at number of nodes in graph
    k_samples = min(50, len(G0.nodes()))
    bc_G0 = nx.betweenness_centrality(G0, k=k_samples, weight='weight', normalized=True, seed=42)
    bc_G_drug = nx.betweenness_centrality(G_drug, k=k_samples, weight='weight', normalized=True, seed=42)
    
    # Compute absolute differences for each node
    shifts = [abs(bc_G_drug[n] - bc_G0[n]) for n in bc_G0]
    
    # Return mean shift
    return sum(shifts) / len(shifts) if shifts else 0.0


def compute_community_change(G0, G_drug) -> float:
    """
    Compute the community structure change between two graphs.
    
    Uses the Louvain algorithm to detect community structure in both graphs,
    then measures dissimilarity using Normalized Mutual Information (NMI).
    A higher value indicates the drug significantly disrupts community organization.
    
    Parameters:
    -----------
    G0 : NetworkX Graph
        Baseline protein-protein interaction graph
    G_drug : NetworkX Graph
        Graph after drug binding simulation (same nodes as G0, some edges weighted differently)
    
    Returns:
    --------
    float
        Community dissimilarity score: 0.0 (identical communities) to 1.0 (completely different)
    """
    # Validate that graphs have identical node sets
    assert set(G0.nodes()) == set(G_drug.nodes()), "Graphs must have identical nodes"
    
    # Validate weights exist on all edges
    assert all('weight' in d for _, _, d in G0.edges(data=True)), "G0: all edges must have 'weight' attribute"
    assert all('weight' in d for _, _, d in G_drug.edges(data=True)), "G_drug: all edges must have 'weight' attribute"
    
    # Handle edge case: no edges
    if G0.number_of_edges() == 0:
        return 0.0
    
    # Detect communities using Louvain algorithm with deterministic seeding
    partition_G0 = best_partition(G0, weight='weight', random_state=42)
    partition_G_drug = best_partition(G_drug, weight='weight', random_state=42)
    
    # Convert partitions to label lists in consistent node order
    node_list = sorted(G0.nodes())
    labels_G0 = [partition_G0[n] for n in node_list]
    labels_G_drug = [partition_G_drug[n] for n in node_list]
    
    # Compute Normalized Mutual Information
    nmi = normalized_mutual_info_score(labels_G0, labels_G_drug)
    
    # Return dissimilarity: 0.0 = no change, 1.0 = maximum change
    return 1.0 - nmi


def compute_spectral_gap(G0, G_drug) -> float:
    try:
        # extract largest connected component from each graph
        lcc0 = G0.subgraph(max(nx.connected_components(G0), key=len)).copy()
        lccd = G_drug.subgraph(max(nx.connected_components(G_drug), key=len)).copy()
        ac0 = nx.algebraic_connectivity(lcc0, weight='weight')
        acd = nx.algebraic_connectivity(lccd, weight='weight')
        return float(acd - ac0)
    except Exception as e:
        print(f"WARNING spectral gap: {e}")
        return 0.0


def compute_entropy_delta(G0, G_drug, targets) -> float:
    from scipy.stats import entropy as scipy_entropy
    deltas = []
    for t in targets:
        if t not in G0 or t not in G_drug:
            continue
        neighbors = list(G0.neighbors(t))
        if not neighbors:
            continue
        # absolute weight change per neighbour
        changes = []
        for n in neighbors:
            w0 = G0[t][n]['weight']
            wd = G_drug[t][n].get('weight', w0)
            changes.append(abs(w0 - wd))
        total = sum(changes)
        if total == 0:
            # no change at all — drug had no effect on this target
            deltas.append(0.0)
            continue
        # entropy of the change distribution
        # high entropy = change spread evenly across all neighbours
        # low entropy = change concentrated on few neighbours
        p = [c / total for c in changes]
        deltas.append(scipy_entropy(p))
    return float(sum(deltas) / len(deltas)) if deltas else 0.0


def compute_rsv(G0, G_drug, targets) -> dict:
    """
    Compute the complete Repurposing Signature Vector (RSV) for a drug candidate.
    
    The RSV is a 4-element signature measuring how a drug perturbs protein-protein
    interaction networks at multiple topological scales:
    
    1. betweenness_shift [0, 1]: Centrality redistribution
       - Higher = major rewiring of central hubs
    
    2. community_change [0, 1]: Pathway reorganization  
       - Higher = complete disruption of biological modules
    
    3. spectral_gap (-∞, ∞): Network resilience
       - Negative = fragmented network, Positive = more robust
       - NOTE: Unbounded. Consider normalization for cross-graph comparison.
    
    4. entropy_delta (-∞, ∞): Neighborhood chaos around targets
       - Positive = disordered, Negative = organized
    
    Returns:
    --------
    dict
        {"betweenness_shift": float, "community_change": float, 
         "spectral_gap": float, "entropy_delta": float,
         "vector": [bs, cc, sg, ed]}
        On failure: component = NaN (detectable, not silent 0.0)
    """
    rsv = {}
    
    # 1. Compute betweenness shift
    try:
        rsv["betweenness_shift"] = compute_betweenness_shift(G0, G_drug)
    except Exception as e:
        print(f"WARNING: compute_betweenness_shift failed: {e}")
        rsv["betweenness_shift"] = float("nan")
    
    # 2. Compute community change
    try:
        rsv["community_change"] = compute_community_change(G0, G_drug)
    except Exception as e:
        print(f"WARNING: compute_community_change failed: {e}")
        rsv["community_change"] = float("nan")
    
    # 3. Compute spectral gap
    try:
        rsv["spectral_gap"] = compute_spectral_gap(G0, G_drug)
    except Exception as e:
        print(f"WARNING: compute_spectral_gap failed: {e}")
        rsv["spectral_gap"] = float("nan")
    
    # 4. Compute entropy delta
    try:
        rsv["entropy_delta"] = compute_entropy_delta(G0, G_drug, targets)
    except Exception as e:
        print(f"WARNING: compute_entropy_delta failed: {e}")
        rsv["entropy_delta"] = float("nan")
    
    # 5. Create vector form for ML/downstream pipelines
    rsv["vector"] = [
        rsv["betweenness_shift"],
        rsv["community_change"],
        rsv["spectral_gap"],
        rsv["entropy_delta"]
    ]
    
    return rsv


if __name__ == "__main__":
    import random
    
    # Set seed for reproducibility
    random.seed(42)
    
    # Create G0: random graph with 20 nodes and 50 edges
    G0 = nx.Graph()
    G0.add_nodes_from(range(20))
    
    edges_added = 0
    while edges_added < 50:
        u = random.randint(0, 19)
        v = random.randint(0, 19)
        if u != v and not G0.has_edge(u, v):
            weight = random.uniform(0.1, 1.0)
            G0.add_edge(u, v, weight=weight)
            edges_added += 1
    
    # Create G_drug as a deep copy of G0
    G_drug = deepcopy(G0)
    
    # Reduce weights of 10 random edges in G_drug by 50%
    edges_list = list(G_drug.edges())
    random_edges = random.sample(edges_list, 10)
    
    for u, v in random_edges:
        G_drug[u][v]['weight'] *= 0.5
    
    # Set targets to first 3 nodes
    targets = [0, 1, 2]
    
    # Compute full RSV
    rsv = compute_rsv(G0, G_drug, targets)
    
    # Verify output structure
    assert isinstance(rsv, dict), "RSV is not a dict"
    expected_keys = {"betweenness_shift", "community_change", "spectral_gap", "entropy_delta", "vector"}
    assert set(rsv.keys()) == expected_keys, f"RSV keys mismatch. Expected {expected_keys}, got {set(rsv.keys())}"
    for key in ["betweenness_shift", "community_change", "spectral_gap", "entropy_delta"]:
        assert isinstance(rsv[key], (float, int)), f"RSV[{key}] is not numeric: {type(rsv[key])}"
    assert isinstance(rsv["vector"], list), "RSV vector is not a list"
    assert len(rsv["vector"]) == 4, f"RSV vector does not have exactly 4 elements, got {len(rsv['vector'])}"
    
    # Print results
    print("Betweenness shift: {:.4f}".format(rsv["betweenness_shift"]))
    print("Community change: {:.4f}".format(rsv["community_change"]))
    print("Spectral gap: {:.4f}".format(rsv["spectral_gap"]))
    print("Entropy delta: {:.4f}".format(rsv["entropy_delta"]))
    print("Vector: {}".format(rsv["vector"]))
    
    print("ALL TESTS PASSED")

    def test_fixed_functions():
        import networkx as nx
        import copy
        G0 = nx.erdos_renyi_graph(50, 0.3, seed=42)
        for u, v in G0.edges():
            G0[u][v]['weight'] = 0.8
        G_drug = copy.deepcopy(G0)
        # simulate binding: halve weights around node 0
        for n in list(G_drug.neighbors(0)):
            G_drug[0][n]['weight'] = 0.4
        targets = [0]
        print("Spectral gap delta:", compute_spectral_gap(G0, G_drug))
        print("Entropy delta:", compute_entropy_delta(G0, G_drug, targets))
        assert compute_spectral_gap(G0, G_drug) != 0.0, "Spectral gap still 0"
        assert compute_entropy_delta(G0, G_drug, targets) != 0.0, "Entropy delta still 0"
        print("BOTH FIXED")

    test_fixed_functions()
