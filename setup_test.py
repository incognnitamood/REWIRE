import networkx as nx
import random
import copy

# Set seed for reproducibility
random.seed(42)

# Create G0: random graph with 20 nodes and 50 edges
G0 = nx.Graph()
G0.add_nodes_from(range(20))

# Add 50 random edges with random weights
edges_added = 0
while edges_added < 50:
    u = random.randint(0, 19)
    v = random.randint(0, 19)
    if u != v and not G0.has_edge(u, v):
        weight = random.uniform(0.1, 1.0)
        G0.add_edge(u, v, weight=weight)
        edges_added += 1

# Create G_drug as a deep copy of G0
G_drug = copy.deepcopy(G0)

# Reduce weights of 10 random edges in G_drug by 50%
edges_list = list(G_drug.edges())
random_edges = random.sample(edges_list, 10)

for u, v in random_edges:
    G_drug[u][v]['weight'] *= 0.5

print(f"Modified edges: {len(random_edges)}")

# Function to print graph stats
def print_stats(graph, name):
    num_nodes = graph.number_of_nodes()
    num_edges = graph.number_of_edges()
    weights = [graph[u][v]['weight'] for u, v in graph.edges()]
    avg_weight = sum(weights) / len(weights) if weights else 0
    print(f"{name}:")
    print(f"  Nodes: {num_nodes}")
    print(f"  Edges: {num_edges}")
    print(f"  Average weight: {avg_weight:.4f}")


# Print stats
print_stats(G0, "G0")
print_stats(G_drug, "G_drug")

# Test independence: modify G_drug and check G0 is unchanged
test_edge = list(G_drug.edges())[0]
u, v = test_edge
original_g0_weight = G0[u][v]['weight']

# Modify G_drug
G_drug[u][v]['weight'] = 0.999

# Check that G0 is unchanged
if G0[u][v]['weight'] == original_g0_weight:
    print(f"\n✓ Graphs are independent")
else:
    print(f"\n✗ Graphs are NOT independent")

print("SETUP OK")
