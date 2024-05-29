import networkx as nx

def max_weight_dfs_tree(G: nx.DiGraph,
                        source: str,
                        reference_path: list = []) -> tuple(nx.DiGraph, nx.DiGraph):
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
    # stack += [(reference_path[i], reference_path[i - 1] if i > 0 else source) for i in range(len(reference_path))]

    for i in range(len(reference_path)):
        if i > 0:
            stack.append((reference_path[i], reference_path[i - 1]))
            ancestors.add(reference_path[i - 1])
        if i > 1:
            dfs_tree.add_edge(reference_path[i - 2], reference_path[i - 1],
                              weight=G[reference_path[i - 2]][reference_path[i - 1]]['weight'])
            ref_dag.add_edge(reference_path[i - 2], reference_path[i - 1],
                             weight=G[reference_path[i - 2]][reference_path[i - 1]]['weight'])

    last_visited_node = reference_path[-2]

    step_counter = 0
    while stack:
        current_node, parent = stack.pop()
        #ic(current_node, parent, last_visited_node, ancestors)
        if current_node not in visited:
            # visit current_node
            G.nodes[current_node]['first_visit'] = step_counter
            step_counter += 1

            # Find path from previous_node to current_node
            if parent != last_visited_node:

                path = nx.shortest_path(dfs_tree, parent, last_visited_node)
                for v in path:
                    if 'second_visit' not in G.nodes[v]:
                        G.nodes[v]['second_visit'] = step_counter
                    G.nodes[v]['last_visit'] = step_counter
                    step_counter += 1
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
            # reference_node = None

            # Last-added edges will be visited first
            for _, neighbor, _ in neighbors:
                # if (current_node, neighbor) in reference_path:
                #     assert(reference_node is None)
                #     reference_node = neighbor
                #     continue
                # (current_node, neighbor) not in stack and
                if current_node != neighbor and neighbor not in ancestors:
                    stack.append((neighbor, current_node))


            # if reference_node is not None:
            #     stack.append((reference_node, current_node))

            last_visited_node = current_node

        elif current_node not in ancestors:
            if parent is not None:
                ref_dag.add_edge(parent, current_node, weight=G[parent][current_node]['weight'])

    return dfs_tree, ref_dag
