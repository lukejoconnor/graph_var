from math import inf
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from utils import flip, read_gfa
from search_tree import max_weight_dfs_tree
from kruskal import kruskal_mst
# from icecream import ic


class DiED_Graph(nx.DiGraph):
    def __init__(self,
                 gfa_file: str = None,
                 edgeinfo_file: str = None,
                 reference_path_index: int = None,
                 spanning_forest_method: str = 'dfs'
                 ):
        super(DiED_Graph, self).__init__()
        self.reference_tree: nx.classes.digraph.DiGraph = nx.DiGraph()
        self.reference_dag: nx.classes.digraph.DiGraph = nx.DiGraph()
        self.variant_edges: set = set()
        self.walks: list = []
        self.reference_path: list = []
        self.position: dict = {}
        self.num_walks: int = 0
        self.num_variants: int = 0
        self.genotypes: list = []

        if gfa_file is None:
            return

        genes, edges, walks = read_gfa(gfa_file)
        print("Num of Genes:", len(genes))
        print("Num of Edges:", len(edges))
        for gene in genes:
            self.add_dinode(gene)
        for n, edge in enumerate(edges):
            self.add_diedge(edge[0], edge[1], edge[2], edge[3], n+1)
            self.add_diedge(edge[1], edge[0], flip(edge[3]), flip(edge[2]), -(n+1))
        self.walks = walks
        self.num_walks = len(walks)

        # Add universal source and sink nodes
        self.add_source_and_sink_nodes()
        print("finished adding source")

        # Add reference path
        # (Done)TODO check that reference path has no duplicate vertices
        if reference_path_index is not None:
            if reference_path_index >= self.num_walks or reference_path_index < 0:
                raise ValueError(f'Reference walk index should be an integer >= 0 and < {self.num_walks}')
            self.reference_path = self.walks[reference_path_index]

        assert len(self.reference_path) == len(set(self.reference_path)), "The reference path has duplicate vertices."

        if edgeinfo_file:
            # Read edgeinfo file and create reference tree and reference dag
            self.read_edgeinfo(edgeinfo_file)

            # Because edges involving start + end nodes are not in the edge_info file, annotate them now
            print("Start annotating")
            self.annotate_start_and_end()

            all_edges = set(self.edges())
            reference_edges = set(self.reference_tree.edges())
            self.variant_edges = all_edges.difference(reference_edges)
        else:
            # Add weights to edges
            print("Start adding weights")
            self.add_weights()
            print("Finish adding weights, start creating tree")
            # Create spanning tree
            self.add_reference(reference_walk_index=reference_path_index, spanning_forest_method=spanning_forest_method)

        # Add positions
        print("Finish creating tree, start adding position")
        self.add_positions()

        # Call variants
        print("Finish adding position, start counting variant")
        self.count_variants()

    def add_source_and_sink_nodes(self): # (Done)TODO add edges between universal source and any start point of a walk (even if not a source); same with universal end
        source_nodes = set([node for node, in_degree in self.in_degree()])
        sink_nodes = set([node for node, out_degree in self.out_degree()])

        self.start_node = 'start_node'
        self.add_node(self.start_node)
        self.end_node = 'end_node'
        self.add_node(self.end_node)

        for source_node in source_nodes:
            self.add_edge(self.start_node, source_node, weight=0)

        for sink_node in sink_nodes:
            self.add_edge(sink_node, self.end_node, weight=0)

        for i in range(len(self.walks)):
            self.walks[i] = ['start_node'] + self.walks[i]
            self.walks[i] = self.walks[i] + ['end_node']
            if not self.has_edge(self.start_node, self.walks[i][1]):
                self.add_edge(self.start_node, self.walks[i][1], weight=1)
            if not self.has_edge(self.walks[i][-2], self.end_node):
                self.add_edge(self.walks[i][-2], self.end_node, weight=1)


    def annotate_start_and_end(self):

        # Add weights to these edges
        for walk in self.walks:
            if len(walk) < 2:
                continue
            first_node = walk[1]
            last_node = walk[-2]
            assert(walk[0] == self.start_node and walk[-1] == self.end_node)

            if first_node in self[self.start_node]:
                self[self.start_node][first_node]['weight'] += 1
            else:
                self.add_edge(self.start_node, first_node, weight=1)

            if self.end_node in self[last_node]:
                self[last_node][self.end_node]['weight'] += 1
            else:
                self.add_edge(last_node, self.end_node, weight=1)

        # All out-edges of start node are in both the ref tree and the ref dag
        for edge in self.out_edges(self.start_node):
            u,v = edge
            self.reference_dag.add_edge(u,v)
            self.reference_tree.add_edge(u,v)

        # All in-edges of the end node are in the ref dag
        for edge in self.in_edges(self.end_node):
            u,v = edge
            self.reference_dag.add_edge(u,v)

        # Exactly one in-edge of the end node is in the ref tree
        if len(self.reference_path) > 1:
            # Pick the reference walk penultimate node
            self.reference_tree.add_edge(self.reference_path[-2], self.end_node)
        else:
            # Pick the max-weight penultimate node
            max_weight_end_edge = max(self.in_edges(self.end_node, data='weight'), key=lambda tup: tup[2])
            self.reference_tree.add_edge(max_weight_end_edge[0], self.end_node)

    def add_reference(self, reference_walk_index: int = None, spanning_forest_method: str = 'dfs'):

        initial_reference_tree = nx.DiGraph()

        # Construct reference tree and ref DAG
        if spanning_forest_method == 'kruskal':
            for u, v in zip(self.reference_path, self.reference_path[1:]):
                # Step 3: Check if the edge exists in G and then add it to H with the same attributes
                if self.has_edge(u, v):
                    edge_attr = self[u][v]  # Get the attributes of the edge from G
                    initial_reference_tree.add_edge(u, v, **edge_attr)  # Add the edge to H with the attributes

            self.reference_tree = kruskal_mst(
                self,
                universal_source=self.start_node,
                initial_graph=initial_reference_tree
            )
        elif spanning_forest_method == 'dfs':
            self.reference_tree, self.reference_dag = max_weight_dfs_tree(self,
                                                                          source=self.start_node,
                                                                          reference_path=self.reference_path)
        else:
            raise ValueError(f"spanning_forest_method should be either 'kruskal' or 'dfs'")

        # All of the edges in the pangenome graph not in the reference tree
        all_edges = set(self.edges())
        reference_edges = set(self.reference_tree.edges())
        self.variant_edges = all_edges.difference(reference_edges)
        self.num_variants = len(self.variant_edges)

    def write_edgeinfo(self, filename: str) -> None:
        with open(filename, 'w') as file:
            # Write the header row
            file.write("index,weight,in_reference_tree,in_reference_dag\n")

            # Sort the edges by 'index' attribute for writing in sorted order
            edges = [edge for edge in self.edges(data=True) if 'index' in edge[2]]
            sorted_edges = sorted(edges, key=lambda edge: edge[2]['index'])

            for edge in sorted_edges:
                u, v, data = edge
                weight = data['weight']
                index = data['index']

                # Check if the edge is in the reference_tree
                in_ref_tree = 1 if self.reference_tree.has_edge(u, v) else 0

                # Check if the edge is in the reference_dag
                in_ref_dag = 1 if self.reference_dag.has_edge(u, v) else 0

                # Write the edge information to the file
                file.write(f"{index},{weight},{in_ref_tree},{in_ref_dag}\n")

    def read_edgeinfo(self, filename: str) -> None:

        weight = np.zeros(2 * self.number_of_edges(), dtype=np.int32)
        in_reference_tree = np.zeros(2 * self.number_of_edges(), dtype=np.int32)
        in_reference_dag = np.zeros(2 * self.number_of_edges(), dtype=np.int32)

        # Read the file
        with open(filename, 'r') as file:
            # Skip the header row
            next(file)

            for line in file:
                # Split the line into components
                parts = line.strip().split(',')
                assert(len(parts) == 4)

                # Extract the edge info
                index = int(parts[0])
                weight[index] = int(parts[1])
                in_reference_tree[index] = int(parts[2])
                in_reference_dag[index] = int(parts[3])

        # Use the edge info to add weights, create ref tree + ref dag
        for edge in self.edges(data=True):
            u, v, data = edge
            if 'index' not in data:
                continue
            idx = self[u][v]['index']
            self[u][v]['weight'] = weight[idx]
            if in_reference_dag[idx] == 0:
                continue # also not in ref tree
            self.reference_dag.add_edge(u, v)
            if in_reference_tree[idx] == 0:
                continue
            self.reference_tree.add_edge(u, v)

    def add_positions(self):
        if self.reference_dag.number_of_nodes() < self.number_of_nodes():
            raise ValueError('Reference DAG does not have enough nodes, probably because it is not defined')

        # Define node positions
        for direction in range(2):
            self.get_node_positions(direction)

        # Check consistency of node positions
        for node in self.reference_dag.nodes():
            assert (self.position[0][node] <= self.position[1][node])

    def add_dinode(self, gene, seq=''):
        self.add_node(str(gene)+'_+', sequence=seq, direction='+')
        self.add_node(str(gene)+'_-', sequence=seq, direction='-')

    def add_diedge(self, gene1, gene2, direction1, direction2, n):
        self.add_edge(str(gene1)+'_'+direction1, str(gene2)+'_'+direction2, weight=0, index=n)

    def add_weights(self):
        for walk_list in self.walks:
            for i in range(int(len(walk_list))-1):
                #print(self.edges)
                self.edges[walk_list[i], walk_list[i+1]]['weight'] += 1

    def count_variants(self):
        self.genotypes = [{} for i in range(self.num_walks)]
        for individual in range(self.num_walks):
            walk_list = self.walks[individual]
            for i in range(len(walk_list) - 1):
                edge = walk_list[i], walk_list[i + 1]
                if edge not in self.variant_edges:
                    continue

                if edge in self.genotypes[individual]:
                    self.genotypes[individual][edge] += 1
                else:
                    self.genotypes[individual][edge] = 1

    def count_edge_visits(self, geno: dict) -> dict:

        sources = {self.start_node: 1}
        sinks = list([self.end_node])

        # For alternative alleles (u,v), add v to the set of "sources" and u to the set of "sinks"
        for variant_edge, visit_count in geno.items():
            u,v = variant_edge
            if v in sources:
                sources[v] += visit_count
            else:
                sources[v] = visit_count # TODO correct?
            for i in range(visit_count):
                sinks.append(u)

        edge_visits = geno.copy()
        for sink in sinks:
            # Walk up the tree (in the only possible direction) until reaching a source
            current_node = sink
            while current_node not in sources:
                # If current_node is the root, it means that the input genotype was invalid
                if current_node == self.start_node:
                    raise ValueError("The input genotype does not correspond to any valid walk")

                # Increment the counter
                previous_node = current_node
                current_node = next(self.reference_tree.predecessors(current_node))
                if (current_node, previous_node) in edge_visits:
                    edge_visits[(current_node, previous_node)] += 1
                else:
                    edge_visits[(current_node, previous_node)] = 1

            # when reaching the source, it is "used up"
            if sources[current_node] == 1:
                sources.pop(current_node)
            else:
                sources[current_node] -= 1

        return edge_visits


    def get_node_positions(self, direction):
        self.position[direction] = {u:-inf * (-1)**direction for u in self.reference_dag.nodes}

        for n, u in enumerate(self.reference_path):
            self.position[direction][u] = n

        if direction == 0:
            G = self.reference_dag
        else:
            G = self.reference_dag.reverse()

        order = nx.topological_sort(G)
        for u in order:
            pred = list(G.predecessors(u))
            pred.append(u)
            if direction == 0:
                self.position[direction][u] = np.max([self.position[direction][v] for v in pred]).astype(int)
            else:
                self.position[direction][u] = np.min([self.position[direction][v] for v in pred]).astype(int)


    def find_path_to_linear_ref(self, bfs_tree, node):
        current = node
        while current not in self.linear_reference + [self.start_node]:
            # In a BFS tree, each node except the root has exactly one predecessor
            preds = list(bfs_tree.predecessors(current))
            current = preds[0]
        return current

    def show(self, G, i_pos=None):
        if i_pos is None:
            pos = nx.spring_layout(G)
        else:
            pos = i_pos

        plt.figure(figsize=(15,10))
        # Draw the graph
        nx.draw_networkx_nodes(G, pos, node_color='white', node_shape='o', edgecolors='black', node_size=500)

        # Draw edges
        for edge, key in enumerate(G.edges):
            if key[2] == 'r':
                nx.draw_networkx_edges(G, pos, arrows=True, arrowstyle='-', edgelist=[key], connectionstyle=f'arc3,rad={0}', edge_color='red')
            else:
                nx.draw_networkx_edges(G, pos, arrows=True, arrowstyle='-', edgelist=[key], connectionstyle=f'arc3,rad={0.30}', edge_color='black')

        # Draw node labels
        nx.draw_networkx_labels(G, pos, font_size=8)

        # Show plot
        plt.show()
