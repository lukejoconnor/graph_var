import networkx as nx

def max_weight_dfs_tree(G, source):
    """
    Perform a depth-first search (DFS) prioritizing higher weight edges and return a DFS tree.

    :param G: A NetworkX directed graph
    :param source: The root node for the DFS tree
    :return: A NetworkX directed graph representing the DFS tree with prioritized higher weights
    """
    visited = set()  # Keep track of visited nodes
    ancestors = set() # Ancestors of last visited node
    last_visted_node = None
    stack = [(source, None)]  # Stack for DFS, storing (node, parent)
    dfs_tree = nx.DiGraph()  # Initialize an empty directed graph for the DFS tree
    ref_dag = nx.DiGraph()

    while stack:
        current_node, parent = stack.pop()
        #ic(current_node, parent, last_visted_node, ancestors)
        if current_node not in visited:
            # visit current_node
            # Find path from previous_node to current_node
            if parent != last_visted_node:
                path = nx.shortest_path(dfs_tree, parent, last_visted_node)
                for v in path:
                    if v is not parent:
                      ancestors.remove(v)

            visited.add(current_node)
            ancestors.add(current_node)
            if parent is not None:
                dfs_tree.add_edge(parent, current_node, weight=G[parent][current_node]['weight'])
                ref_dag.add_edge(parent, current_node, weight=G[parent][current_node]['weight'])
            # Sort neighbors by descending weight
            neighbors = sorted(G.edges(current_node, data=True), key=lambda x: x[2]['weight'], reverse=False)
            for _, neighbor, _ in neighbors:
                stack.append((neighbor, current_node))
            last_visted_node = current_node

        elif current_node not in ancestors:
            if parent is not None:
                ref_dag.add_edge(parent, current_node, weight=G[parent][current_node]['weight'])



    return dfs_tree, ref_dag