import networkx as nx

def find(parent, i):
    if parent[i] == i:
        return i
    return find(parent, parent[i])

def union(parent, rank, x, y):
    xroot = find(parent, x)
    yroot = find(parent, y)

    if rank[xroot] < rank[yroot]:
        parent[xroot] = yroot
    elif rank[xroot] > rank[yroot]:
        parent[yroot] = xroot
    else:
        parent[yroot] = xroot
        rank[xroot] += 1

def kruskal_mst(graph, universal_source = None, initial_graph = nx.DiGraph()):
    result = []  # This will store the resultant MST
    i, e = 0, 0  # Number of edges to be taken is equal to V-1

    # Step 1: Sort all the edges in non-decreasing order of their weight
    initial_edges = set(initial_graph.edges())
    edges = [edge for edge in graph.edges(data='weight') if not edge in initial_edges]
    edges = sorted(edges, key=lambda t: t[2], reverse=True)

    F = initial_graph.copy()
    F.add_nodes_from(graph.nodes())

    parent, rank, in_degree = {}, {}, {}
    for node in graph.nodes():
        parent[node] = node
        rank[node] = 0

    # Number of edges to be taken is equal to V-1
    for i in range(len(edges)):
        # Step 2: Pick the smallest edge. Check if it forms a cycle with the spanning tree formed so far.
        u, v, w = edges[i]

        if F.in_degree(v) > 0:
            continue

        x = find(parent, u)
        y = find(parent, v)

        # If including this edge does not cause a cycle, include it in the result and increment the index of the result for the next edge.
        if x != y:
            e = e + 1
            F.add_edge(u, v, weight=w)
            union(parent, rank, x, y)

    if universal_source is not None:
        source_nodes = [node for node, in_degree in F.in_degree() if in_degree == 0]
        for source_node in source_nodes:
            if source_node == universal_source:
                continue
            F.add_edge(universal_source, source_node, weight=0)

    return F

# Example usage:
if __name__ == '__main__':
    G = nx.DiGraph()
    G.add_edge(2, 3, weight=4)
    G.add_edge(3, 4, weight=3)
    G.add_edge(4, 1, weight=4)
    G.add_edge(1, 3, weight=5)

    F = kruskal_mst(G)

    print(F.edges())