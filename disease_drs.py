import csv
import json

import networkx as nx
import numpy as np
import requests
from sklearn.metrics.cluster import normalized_mutual_info_score


def get_disease_genes(disease_name: str) -> list:
    query = """
    query DiseaseGenes($disease: String!) {
      search(queryString: $disease, entityNames: ["disease"], page: {index: 0, size: 1}) {
        hits {
          id
          name
        }
      }
      disease(id: $diseaseId) {
        associatedTargets(page: {index: 0, size: 50}) {
          rows {
            target {
              id
              approvedSymbol
            }
            score
          }
        }
      }
    }
    """

    search_query = """
    query DiseaseSearch($disease: String!) {
      search(queryString: $disease, entityNames: ["disease"], page: {index: 0, size: 1}) {
        hits {
          id
          name
        }
      }
    }
    """

    search_payload = {"query": search_query, "variables": {"disease": disease_name}}
    search_resp = requests.post(
        "https://api.platform.opentargets.org/api/v4/graphql",
        headers={"Content-Type": "application/json"},
        data=json.dumps(search_payload),
        timeout=30,
    )
    search_resp.raise_for_status()
    search_data = search_resp.json()
    hits = search_data.get("data", {}).get("search", {}).get("hits", [])
    if not hits:
        print(f"{disease_name}: 0 genes returned")
        return []

    disease_id = hits[0]["id"]

    main_payload = {
        "query": query,
        "variables": {"disease": disease_name, "diseaseId": disease_id},
    }
    main_resp = requests.post(
        "https://api.platform.opentargets.org/api/v4/graphql",
        headers={"Content-Type": "application/json"},
        data=json.dumps(main_payload),
        timeout=30,
    )
    main_resp.raise_for_status()
    main_data = main_resp.json()
    rows = (
        main_data.get("data", {})
        .get("disease", {})
        .get("associatedTargets", {})
        .get("rows", [])
    )

    genes = [
        {"symbol": row["target"]["approvedSymbol"], "score": row["score"]}
        for row in rows
        if row.get("target") and row["target"].get("approvedSymbol")
    ]

    print(f"{disease_name}: {len(genes)} genes returned")
    return genes


  def load_ppi_graph(ppi_path: str) -> nx.Graph:
    graph = nx.Graph()
    with open(ppi_path, "r", encoding="utf-8") as handle:
      reader = csv.DictReader(handle)
      for row in reader:
        node_a = row.get("gene1") or row.get("protein1") or row.get("node1")
        node_b = row.get("gene2") or row.get("protein2") or row.get("node2")
        weight_raw = row.get("weight") or row.get("combined_score")
        if not node_a or not node_b or weight_raw is None:
          continue
        try:
          weight = float(weight_raw)
        except ValueError:
          continue
        graph.add_edge(node_a, node_b, weight=weight)
    return graph


  def shannon_entropy(values: list) -> float:
    if not values:
      return 0.0
    total = float(sum(values))
    if total <= 0:
      return 0.0
    probs = [v / total for v in values if v > 0]
    if not probs:
      return 0.0
    return float(-sum(p * np.log2(p) for p in probs))


  def community_labels(graph: nx.Graph) -> dict:
    communities = nx.algorithms.community.louvain_communities(
      graph, weight="weight", seed=42
    )
    labels = {}
    for idx, members in enumerate(communities):
      for node in members:
        labels[node] = idx
    return labels


  def algebraic_connectivity(graph: nx.Graph) -> float:
    if graph.number_of_nodes() < 2:
      return 0.0
    laplacian = nx.laplacian_matrix(graph, weight="weight").toarray()
    eigenvalues = np.linalg.eigvalsh(laplacian)
    if eigenvalues.size < 2:
      return 0.0
    return float(np.sort(eigenvalues)[1])


  def compute_drs(disease_name: str) -> np.ndarray:
    genes = get_disease_genes(disease_name)
    if not genes:
      print(f"WARNING: No genes found for {disease_name}; using zero vector")
      return np.zeros(4, dtype=float)

    try:
      base_graph = load_ppi_graph("ppi.csv")
    except FileNotFoundError:
      print("ERROR: ppi.csv not found; cannot compute DRS")
      return np.zeros(4, dtype=float)

    disease_graph = base_graph.copy()
    gene_scores = {entry["symbol"]: float(entry["score"]) for entry in genes}

    for gene, score in gene_scores.items():
      if gene not in disease_graph:
        continue
      for neighbor in list(disease_graph.neighbors(gene)):
        edge_data = disease_graph[gene][neighbor]
        weight = float(edge_data.get("weight", 0.0))
        edge_data["weight"] = weight * (1.0 - score)

    base_bc = nx.betweenness_centrality(base_graph, weight="weight")
    disease_bc = nx.betweenness_centrality(disease_graph, weight="weight")
    bc_changes = [
      abs(disease_bc[node] - base_bc[node]) for node in base_graph.nodes()
    ]
    metric_0 = float(np.mean(bc_changes)) if bc_changes else 0.0

    base_labels = community_labels(base_graph)
    disease_labels = community_labels(disease_graph)
    nodes = list(base_graph.nodes())
    base_vec = [base_labels.get(node, -1) for node in nodes]
    disease_vec = [disease_labels.get(node, -1) for node in nodes]
    metric_1 = 1.0 - float(normalized_mutual_info_score(base_vec, disease_vec))

    metric_2 = algebraic_connectivity(disease_graph) - algebraic_connectivity(
      base_graph
    )

    entropy_changes = []
    for gene in gene_scores:
      if gene not in base_graph:
        continue
      base_weights = [
        float(base_graph[gene][nbr].get("weight", 0.0))
        for nbr in base_graph.neighbors(gene)
      ]
      disease_weights = [
        float(disease_graph[gene][nbr].get("weight", 0.0))
        for nbr in disease_graph.neighbors(gene)
      ]
      entropy_changes.append(
        shannon_entropy(disease_weights) - shannon_entropy(base_weights)
      )
    metric_3 = float(np.mean(entropy_changes)) if entropy_changes else 0.0

    return np.array([metric_0, metric_1, metric_2, metric_3], dtype=float)


if __name__ == "__main__":
  diseases = [
    "Parkinson's Disease",
    "Alzheimer's Disease",
    "Type 2 Diabetes",
    "Breast Cancer",
    "Hypertension",
  ]

  drs_vectors = []
  for disease in diseases:
    drs_vector = compute_drs(disease)
    drs_vectors.append(drs_vector)
    formatted = ", ".join(f"{value:.4f}" for value in drs_vector)
    print(f"{disease}: [{formatted}]")

  drs_matrix = np.vstack(drs_vectors)
  np.save("disease_drs.npy", drs_matrix)
  with open("disease_names.txt", "w", encoding="utf-8") as handle:
    handle.write("\n".join(diseases))

  print("Saved disease_drs.npy shape (5, 4) and disease_names.txt")

  loaded = np.load("disease_drs.npy")
  print(f"Loaded disease_drs.npy shape {loaded.shape}")
  for disease, row in zip(diseases, loaded):
    formatted = ", ".join(f"{value:.4f}" for value in row)
    print(f"{disease}: [{formatted}]")
    if np.allclose(row, 0.0):
      print(f"WARNING: zero vector for {disease}")

  distinct = True
  for i in range(len(loaded)):
    for j in range(i + 1, len(loaded)):
      if np.allclose(loaded[i], loaded[j]):
        distinct = False

  if distinct and not any(np.allclose(row, 0.0) for row in loaded):
    print("Validation passed — all 5 DRS vectors are non-zero and distinct")
