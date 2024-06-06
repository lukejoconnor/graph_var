import networkx as nx


def max_weight_dfs_tree(G: nx.DiGraph,
                        source: str,
                        reference_path: list) -> nx.DiGraph:
    """
    Perform a depth-first search (DFS) prioritizing higher weight edges and return a DFS tree.

    :param G: A NetworkX directed graph
    :param source: The root node for the DFS tree
    :param reference_path: A list of nodes in the reference path
    :return: A NetworkX directed graph representing the DFS tree with prioritized higher weights
    """


    dfs_tree = nx.DiGraph()
    visited = set()  # Keep track of visited nodes
    stack = []  # Stack for DFS, storing (node, parent)

    def visit_node(parent, current_node):
        # visit current_node
        visited.add(current_node)
        if parent is not None:
            dfs_tree.add_edge(parent, current_node)

        # Sort neighbors by descending edge weight
        neighbors = sorted(G.edges(current_node, data=True), key=lambda x: x[2]['weight'], reverse=False)

        # Last-added edges will be visited first
        for _, neighbor, _ in neighbors:
            if neighbor not in visited:
                stack.append((neighbor, current_node))

    for i in range(len(reference_path)):
            if i == 0:
                reference_edge = (None, reference_path[0])
            else:
                reference_edge = (reference_path[i - 1], reference_path[i])
            visit_node(*reference_edge)

    while stack:
        current_node, parent = stack.pop()
        if current_node not in visited:
            visit_node(parent, current_node)

    return dfs_tree


