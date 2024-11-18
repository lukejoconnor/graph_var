import networkx as nx
from utils import node_complement

def max_weight_dfs_tree(G: nx.DiGraph,
                        reference_path: list) -> nx.DiGraph:
    """
    Perform a depth-first search (DFS) prioritizing higher weight edges and return a DFS tree.
    :param G: A NetworkX directed graph
    :param reference_path: A list of nodes in the reference path, beginning with source node
    :return: A DFS spanning tree of G
    """
    dfs_tree = nx.DiGraph()
    visited = set()
    stack = []
    def visit_node(parent, current_node):
        # visit current_node
        visited.add(current_node)
        if parent is not None:
            assert G.has_edge(parent, current_node)
            dfs_tree.add_edge(parent, current_node)

        # Sort neighbors by descending edge weight
        neighbors = sorted(G.edges(current_node, data=True), key=lambda x: x[2]['weight'], reverse=False)

        # Last-added edges will be visited first
        for _, neighbor, _ in neighbors:
            if neighbor not in visited:
                stack.append((current_node, neighbor))

    for i in range(len(reference_path)):
        if i == 0:
            reference_edge = (None, reference_path[0])
        else:
            reference_edge = (reference_path[i - 1], reference_path[i])
        visit_node(*reference_edge)

    while stack:
        parent, current_node = stack.pop()
        if current_node not in visited:
            visit_node(parent, current_node)

    return dfs_tree

def assign_node_directions(G: nx.DiGraph,
                        reference_path: list) -> None:
    """
    Perform a depth-first search (DFS) prioritizing higher weight edges and return a DFS tree.

    :param G: A NetworkX directed graph
    :param reference_path: A list of nodes in the reference path, beginning with source node
    :return: A DFS spanning tree of G
    """

    G_undirected = G.to_undirected(as_view=True)
    visited = set()
    stack = []

    def visit_node(parent, current_node, current_direction):

        # visit current_node
        visited.add(current_node)
        G.nodes[current_node]['direction'] = current_direction
        if parent is not None:
            assert G_undirected.has_edge(parent, current_node)

        # Sort neighbors by descending edge weight
        neighbors = sorted(G_undirected.edges(current_node, data=True), key=lambda x: x[2]['weight'], reverse=False)

        # Last-added edges will be visited first
        for u, neighbor, _ in neighbors:
            if neighbor in visited or node_complement(neighbor) in visited:
                continue
            stack.append((current_node, neighbor, current_direction))
            stack.append((None, node_complement(neighbor), -current_direction))

    for i in range(len(reference_path)):
        if i == 0:
            reference_edge = (None, reference_path[0])
        else:
            reference_edge = (reference_path[i - 1], reference_path[i])
        visit_node(*reference_edge, 1)

    while stack:
        parent, current_node, current_direction = stack.pop()
        if current_node not in visited:
            visit_node(parent, current_node, current_direction)

