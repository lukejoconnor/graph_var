from math import inf
import networkx as nx
import numpy as np
from utils import _flip, read_gfa, _node_complement, _edge_complement
from search_tree import max_weight_dfs_tree
import os


class PangenomeGraph(nx.DiGraph):
    reference_tree: nx.classes.digraph.DiGraph
    reference_path: list[str]
    variant_edges: set
    number_of_biedges: int  # Needed because nx.number_of_edges() runs in O(number of edges)

    # A valid walk proceeds from either +_terminus_+ to -_terminus_+ or from -_terminus_- to +_terminus_-
    @property
    def termini(self) -> tuple[str, str]:
        return '+_terminus', '-_terminus'

    @property
    def biedge_representatives(self) -> list[str]:
        edges = [edge_with_data for edge_with_data in self.edges(data=True) if edge_with_data[2]['is_representative']]
        return sorted(edges, key=lambda edge: edge[2]['index'])

    @property
    def biedge_attribute_names(self) -> tuple:
        return 'index', 'weight', 'is_in_tree', 'is_back_edge'
        # TODO add 'is_in_reference_path' such that we can load with edgeinfo file not specifying reference walk index

    @property
    def node_attribute_names(self) -> tuple:
        return 'direction', 'sequence', 'position'

    def __init__(self,
                 directed_graph: nx.classes.digraph.DiGraph = nx.DiGraph(),
                 reference_tree: nx.classes.digraph.DiGraph = nx.DiGraph(),
                 reference_path: list[str] = None,
                 variant_edges: set = None
                 ):

        super().__init__(directed_graph)
        self.reference_tree = reference_tree
        self.reference_path = reference_path if reference_path else []
        self.variant_edges = variant_edges if variant_edges else {}
        self.number_of_biedges = np.sum(
            [count_or_not for _, _, count_or_not in directed_graph.edges(data='is_representative')]
        )

    @classmethod
    def from_gfa(cls,
                 gfa_file: str,
                 edgeinfo_file: str = None,
                 reference_path_index: int = None,
                 return_walks: bool = False
                 ):

        if not os.path.exists(gfa_file):
            raise FileNotFoundError(gfa_file)
        if edgeinfo_file:
            if not os.path.exists(edgeinfo_file):
                raise FileNotFoundError(edgeinfo_file)

        nodes, edges, walks, sequences = read_gfa(gfa_file)

        print("Num of Nodes:", len(nodes))
        print("Num of Edges:", len(edges))

        # Initialize an instance of the class
        G = cls()

        for node, sequence in zip(nodes, sequences):
            G.add_binode(node, sequence)

        for edge in edges:
            G.add_biedge(edge[0] + '_' + edge[2], edge[1] + '_' + edge[3])

        # Add universal source and sink nodes, also adding terminal edges for each walk
        walk_start_nodes = [walk[0] for walk in walks]
        walk_end_nodes = [walk[-1] for walk in walks]
        G.add_terminal_nodes(walk_start_nodes=walk_start_nodes, walk_end_nodes=walk_end_nodes)

        print("Finished adding start/end nodes")

        # Add reference path
        if reference_path_index is not None:
            if reference_path_index >= len(walks) or reference_path_index < 0:
                raise ValueError(f'Reference walk index should be an integer >= 0 and < {G.num_walks}')

            G.add_reference_path([G.termini[0]+'_+'] + walks[reference_path_index] + [G.termini[1]+'_+'])

        else:
            G.reference_path = nx.shortest_path(G, source=G.termini[0]+'_+', target=G.termini[1]+'_+')

        assert len(G.reference_path) == len(set(G.reference_path)), "The reference path has duplicate vertices."

        if edgeinfo_file:
            # Read edgeinfo file and create reference tree and reference dag
            G.read_edgeinfo(edgeinfo_file)
        else:
            # Add weights to edges
            print("Start adding weights")
            G.compute_edge_weights(walks)
            print("Finish adding weights, start creating tree")

            # Create spanning tree and define variant edges
            G.compute_reference_tree()

        # Add positions
        print("Finish creating tree, start adding position")
        G.compute_node_positions()

        if return_walks:
            return G, walks

        return G

    def add_reference_path(self, reference_path: list):
        self.reference_path = reference_path

        # Define direction of the reference path to be positive
        for node in reference_path:
            if self.nodes[node]['direction'] == -1:
                self.nodes[node]['direction'] = 1
                self.nodes[_node_complement(node)]['direction'] = -1

        # Ensure first and last reference path edges are in the graph
        if not self.has_edge(reference_path[0], reference_path[1]):
            self.add_biedge(reference_path[0], reference_path[1])

        if not self.has_edge(reference_path[-2], reference_path[-1]):
            self.add_biedge(reference_path[-2], reference_path[-1])

    # TODO update find_snps
    # TODO is there a definition of ref/alt that yields a SNP annotation immediately?
    def find_snps(self):
        variant_length = {v: self.position[1][v[1]] - self.position[0][v[0]] for v in self.variant_edges}
        snps = [key for key, value in variant_length.items() if value == 2]
        alt_nodes = [var_edge[0] for var_edge in snps]
        ref_nodes = [self.reference_path[self.position[0][alt_node] + 1] for alt_node in alt_nodes]

        valid_snps = {}
        for idx, (alt_node, ref_node) in enumerate(zip(alt_nodes, ref_nodes)):
            if alt_node != 'start_node':
                alt_seq = self.sequences[alt_node.split('_')[0]]
                ref_seq = self.sequences[ref_node.split('_')[0]]
                if len(alt_seq) == 1 and len(ref_seq) == 1:
                    ref_pos = self.position[0][alt_node] + 1
                    valid_snps[snps[idx]] = (ref_pos, ref_seq, alt_seq)

        return valid_snps

    # TODO annotate variant edges as either dup, del, replacement depending on DAG (may need to modify DFS algorithm
    #  to distinguish forward from crossing edges)
    # TODO add ref/alt as a dict mapping from each SNP edge to a tuple
    def annotate_variants(self):
        return


    def add_terminal_nodes(self, walk_start_nodes: list[str]=None, walk_end_nodes: list[str]=None):

        source_nodes = {node for node, degree in self.in_degree()
                            if degree == 0 and self.nodes[node]['direction'] == 1}

        sink_nodes = {node for node, degree in self.out_degree()
                          if degree == 0 and self.nodes[node]['direction'] == 1}

        if walk_start_nodes:
            source_nodes = source_nodes.union(
                {node for node in walk_start_nodes if self.nodes[node]['direction'] == 1}
            )
            sink_nodes = sink_nodes.union(
                {_node_complement(node) for node in walk_start_nodes if self.nodes[node]['direction'] == -1}
            )
        if walk_end_nodes:
            sink_nodes = sink_nodes.union(
                {node for node in walk_end_nodes if self.nodes[node]['direction'] == 1}
            )
            source_nodes = source_nodes.union(
                {_node_complement(node) for node in walk_end_nodes if self.nodes[node]['direction'] == -1}
            )

        start_node, end_node = self.termini
        self.add_binode(start_node)
        self.add_binode(end_node)

        for source_node in source_nodes:
            self.add_biedge(start_node + '_+', source_node)

        for sink_node in sink_nodes:
            self.add_biedge(sink_node, end_node + '_+')

    def compute_reference_tree(self):

        # Construct reference tree and ref DAG
        positive_subgraph = self.subgraph([n for n, attr in self.nodes(data=True) if attr.get('direction') == 1])
        self.reference_tree, dfs_dag = max_weight_dfs_tree(positive_subgraph,
                                                                      source=self.termini[0] + '_+',
                                                                      reference_path=self.reference_path)

        # Annotate edges
        for edge in positive_subgraph.edges():
            self.edges[edge]['is_in_tree'] = self.reference_tree.has_edge(*edge)
            self.edges[_edge_complement(edge)]['is_in_tree'] = self.reference_tree.has_edge(*edge)
            self.edges[edge]['is_back_edge'] = not dfs_dag.has_edge(*edge)
            self.edges[_edge_complement(edge)]['is_back_edge'] = not dfs_dag.has_edge(*edge)

        # Define variant edge representatives
        self.variant_edges = {(u, v) for u, v, data in self.edges(data=True)
                if data['is_representative'] and not data['is_in_tree']}

    def write_edgeinfo(self, filename: str) -> None:
        with open(filename, 'w') as file:
            # Write the header row
            file.write(','.join(self.biedge_attribute_names) + '\n')

            # Sort the edges by 'index' attribute, restricting to representatives
            sorted_edges = self.biedge_representatives

            for n, edge in enumerate(sorted_edges):
                _, _, data = edge
                assert n == data['index'], 'Something is wrong with edge indices'
                # Write the edge information to the file; map True->'1', False->'0'
                edge_data_list = [str(int(data[key])) for key in self.biedge_attribute_names]
                file.write(','.join(edge_data_list) + '\n')

    def read_edgeinfo(self, filename: str) -> None:

        sorted_edges = self.biedge_representatives

        # Read the file
        with open(filename, 'r') as file:
            header = next(file).strip().split(',')
            for n, line in enumerate(file):
                # Split the line into components
                parts = line.strip().split(',')

                # Find the edge represenative with index == n
                edge = sorted_edges[n][:-1]
                assert self.edges[edge]['index'] == n

                # Extract the edge info
                for i, key in enumerate(self.biedge_attribute_names):
                    self.edges[edge][key] = int(parts[i])
                    self.edges[_edge_complement(edge)][key] = int(parts[i])

        # Define variant edge representatives
        self.variant_edges = {(u, v) for u, v, data in self.edges(data=True)
                              if data['is_representative'] and not data['is_in_tree']}

        # Define reference tree
        for edge in self.edges(data='is_in_tree'):
            if not edge[2]:  # not in tree
                continue

            # Only include forward direction nodes
            if self.nodes[edge[0]]['direction'] == 1:
                self.reference_tree.add_edge(*edge[:-1])

    def add_binode(self, node: str, seq: str = ''):
        node_data = {key: 0 for key in self.node_attribute_names}
        node_data['sequence'] = seq
        node_data['direction'] = 1

        self.add_node(str(node)+'_+', **node_data)

        node_data['sequence'] = reversed(seq)# TODO map bases to complements
        node_data['direction'] = -1

        self.add_node(str(node)+'_-', **node_data)

    def add_biedge(self, node1: str, node2: str, *, weight: int = 0):

        if self.has_edge(node1, node2):
            raise ValueError(f'Attempted to add duplicate biedge: {node1}, {node2}')

        edge_data = {key: 0 for key in self.biedge_attribute_names}
        edge_data['weight'] = weight
        edge_data['index'] = self.number_of_biedges
        edge_data['is_representative'] = False
        self.add_edge(_node_complement(node2), _node_complement(node1), **edge_data)
        edge_data['is_representative'] = True
        self.add_edge(node1, node2, **edge_data)
        self.number_of_biedges += 1

    def representative_edge(self, edge: tuple):
        return edge if self.edges[edge]['is_representative'] else _edge_complement(edge)

    def compute_edge_weights(self, walks):
        for walk in walks:

            for u, v in zip(walk[:-1], walk[1:]):
                self.edges[u, v]['weight'] += 1
                self.edges[_node_complement(v), _node_complement(u)]['weight'] += 1

    def genotype(self, walk: list[str]) -> dict:

        # Append start and end nodes to walk
        start = [self.termini[0] + '_+' if self.nodes[walk[0]]['direction'] == 1 else self.termini[1] + '_-']
        end = [self.termini[1] + '_+' if self.nodes[walk[-1]]['direction'] == 1 else self.termini[0] + '_-']
        walk = start + walk + end

        genotype = {}
        for e in zip(walk[:-1], walk[1:]):
            if not self.has_edge(*e):
                raise ValueError(f"Specified list contains edge {e} which is not present in the graph")

            if self.edges[e]['is_in_tree']:
                continue

            if not self.edges[e]['is_representative']:
                e = _edge_complement(e)

            if e in genotype:
                genotype[e] += 1
            else:
                genotype[e] = 1

        return genotype

    def count_edge_visits(self, genotype: dict) -> dict:

        sinks: list[str] = []

        # Where walk starts will depend on the parity of the number of inversions
        sources: dict[str, int] = {}

        # Add variant edge endpoints as sources or sinks depending on their respective directions
        for variant_edge, visit_count in genotype.items():

            if variant_edge not in self.variant_edges:
                raise ValueError("geno dictionary contains a key which is not a variant edge")

            u, v = variant_edge
            dir_u, dir_v = self.nodes[u]['direction'], self.nodes[v]['direction']

            new_sinks = []
            new_sources = []
            if dir_u == 1:
                new_sinks.append(u)
            else:
                new_sources.append(_node_complement(u))
            if dir_v == 1:
                new_sources.append(v)
            else:
                new_sinks.append(_node_complement(v))

            for w in new_sources:
                if w in sources:
                    sources[w] += visit_count
                else:
                    sources[w] = visit_count

            sinks += visit_count * new_sinks

        # Add sink and source nodes for the beginning and end of the walk, depending on the number of inversions
        num_sources_minus_sinks = np.sum([val for _, val in sources.items()]) - len(sinks)

        # Equal number of + to - and - to + inversions: walk from + terminus to - terminus
        if num_sources_minus_sinks == 0:
            sources[self.termini[0] + '_+'] = 1
            sinks += [self.termini[1] + '_+']

        # Odd number of inversions, with one more from + to - strand: walk from - to -
        elif num_sources_minus_sinks == 2:
            sinks += 2 * [self.termini[1] + '_+']

        # Odd number of inversions, with one more from - to + strand: walk from + to +
        elif num_sources_minus_sinks == -2:
            sources[self.termini[0] + '_+'] = 2

        else:
            raise ValueError("The input genotype does not correspond to any valid walk")


        edge_visits = genotype.copy()
        for sink in sinks:
            # Walk up the tree (in the only possible direction) until reaching a source
            current_node = sink
            while current_node not in sources:
                # If current_node is the root, it means that the input genotype was invalid
                if self.reference_tree.in_degree(current_node) == 0:
                    raise ValueError("The input genotype does not correspond to any valid walk")

                previous_node = current_node
                current_node = next(self.reference_tree.predecessors(current_node))
                edge_representative = self.representative_edge((current_node, previous_node))
                if edge_representative in edge_visits:
                    edge_visits[edge_representative] += 1
                else:
                    edge_visits[edge_representative] = 1

            # when reaching the source, it is "used up"
            if sources[current_node] == 1:
                sources.pop(current_node)
            else:
                sources[current_node] -= 1

        return edge_visits

    def compute_node_positions(self):
        # TODO should position be w.r.t. node count or sequence length?
        # TODO also compute forward position? Needed for anything except defining SNPs?
        for node in self.reference_tree.nodes():
            self.nodes[node]['position'] = -inf

        for n, u in enumerate(self.reference_path):
            self.nodes[u]['position'] = n

        order = nx.topological_sort(self.reference_tree)
        next(order)  # skip the root
        for u in order:
            predecessor = next(self.reference_tree.predecessors(u))
            self.nodes[u]['position'] = np.maximum(self.nodes[u]['position'],
                                               self.nodes[predecessor]['position'])
            self.nodes[_node_complement(u)]['position'] = self.nodes[u]['position']
