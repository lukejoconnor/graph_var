import networkx as nx

def max_weight_dfs_tree(G: nx.DiGraph,
                        source: str,
                        reference_path: list = []) -> nx.DiGraph:
    """
    Perform a depth-first search (DFS) prioritizing higher weight edges and return a DFS tree.

    :param G: A NetworkX directed graph
    :param source: The root node for the DFS tree
    :return: A NetworkX directed graph representing the DFS tree with prioritized higher weights
    """

    dfs_tree = nx.DiGraph()
    ref_dag = nx.DiGraph()

    visited = set()  # Keep track of visited nodes
    ancestors = set()  # Ancestors of last visited node
    last_visited_node = None

    # TODO add entire reference path, in order; then delete special handling
    stack = [(source, None)]  # Stack for DFS, storing (node, parent)

    while stack:
        current_node, parent = stack.pop()
        #ic(current_node, parent, last_visited_node, ancestors)
        if current_node not in visited:
            # visit current_node
            # Find path from previous_node to current_node
            if parent != last_visited_node:
                path = nx.shortest_path(dfs_tree, parent, last_visited_node)
                for v in path:
                    if v is not parent:
                      ancestors.remove(v)

            visited.add(current_node)
            ancestors.add(current_node)
            if parent is not None:
                dfs_tree.add_edge(parent, current_node, weight=G[parent][current_node]['weight'])
                ref_dag.add_edge(parent, current_node, weight=G[parent][current_node]['weight'])

            # Sort neighbors by descending edge weight
            neighbors = sorted(G.edges(current_node, data=True), key=lambda x: x[2]['weight'], reverse=False)

            # Special handling for any edges contained in the reference path
            reference_node = None

            # Last-added edges will be visited first
            for _, neighbor, _ in neighbors:
                if (current_node, neighbor) in reference_path:
                    assert(reference_node is None)
                    reference_node = neighbor
                    continue
                stack.append((neighbor, current_node))

            if reference_node is not None:
                stack.append((reference_node, current_node))

            last_visited_node = current_node

        elif current_node not in ancestors:
            if parent is not None:
                ref_dag.add_edge(parent, current_node, weight=G[parent][current_node]['weight'])



    return dfs_tree, ref_dag