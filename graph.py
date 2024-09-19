from math import inf
import networkx as nx
import numpy as np
from utils import read_gfa, _node_complement, _edge_complement, _sequence_reversed_complement, _node_recover
from search_tree import assign_node_directions, max_weight_dfs_tree
import os
from collections import defaultdict, Counter
from tqdm import tqdm


class PangenomeGraph(nx.DiGraph):
    reference_tree: nx.classes.digraph.DiGraph
    reference_path: list[str]
    variant_edges: set
    number_of_biedges: int  # Needed because nx.number_of_edges() runs in O(number of edges)

    # A valid walk proceeds from either +_terminus_+ to -_terminus_+ or from -_terminus_- to +_terminus_-
    @property
    def termini(self) -> tuple[str, str]:
        return '+_terminus', '-_terminus'

    def is_terminal_node(self, node: str) -> bool:
        return node[:-2] in self.termini

    # Each biedge has a representative edge, whichever is in the .gfa file
    @property
    def sorted_biedge_representatives(self) -> list[str]:
        edges = [edge_with_data for edge_with_data in self.edges(data=True) if edge_with_data[2]['is_representative']]
        return sorted(edges, key=lambda edge: edge[2]['index'])

    @property
    def sorted_positive_variant_edge(self) -> list[str]:
        edges = [edge if self.nodes[edge[0]]['direction'] == 1 or self.is_inversion(edge) else _edge_complement(edge)
                 for edge in list(self.variant_edges)]
        return sorted(edges)

    @property
    def sorted_biedge(self) -> list[str]:
        edges = [edge_with_data for edge_with_data in self.edges(data=True)]
        return sorted(edges, key=lambda edge: edge[2]['index'])

    @property
    def biedge_attribute_names(self) -> tuple:
        return 'index', 'weight', 'is_in_tree', 'branch_point'

    @property
    def node_attribute_names(self) -> tuple:
        return 'direction', 'sequence', 'position', 'forward_position', 'on_reference_path'

    @property
    def vcf_attribute_names(self) -> tuple:
        return 'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'

    @property
    def number_of_binodes(self) -> int:
        return self.number_of_nodes() // 2

    def edge_on_reference_path(self, edge: tuple):
        if not self.has_edge(edge):
            raise ValueError("Graph does not have edge {edge}")
        return self.nodes[edge[1]]['on_reference_path']

    @classmethod
    def from_gfa(cls,
                 gfa_file: str,
                 reference_path_index: int = None,
                 edgeinfo_file: str = None,
                 nodeinfo_file: str = None,
                 return_walks: bool = False,
                 compressed: bool = False
                 ):
        """
        Reads a .gfa file into a PangenomeGraph object.
        """

        if not os.path.exists(gfa_file):
            raise FileNotFoundError(gfa_file)
        if edgeinfo_file:
            if not os.path.exists(edgeinfo_file):
                raise FileNotFoundError(edgeinfo_file)

        data_dict = read_gfa(gfa_file, compressed=compressed)

        binodes = data_dict['nodes']
        biedges = data_dict['edges']
        walks = data_dict['walks']
        walk_sample_names = data_dict['walk_sample_names']
        sequences = data_dict['sequences']
        reference_index = data_dict['reference_index']

        if reference_path_index is not None:
            reference_index = reference_path_index

        print("Reference walk:", reference_index)
        print("Num of binodes:", len(binodes))
        print("Num of biedges:", len(biedges))

        # Initialize an instance of the class
        G = cls()

        for binode, sequence in zip(binodes, sequences):
            G.add_binode(binode, sequence)

        for biedge in biedges:
            node1 = biedge[0] + '_' + biedge[2]
            node2 = biedge[1] + '_' + biedge[3]
            G.add_biedge(node1, node2)

        if reference_index >= len(walks) or reference_index < 0:
            raise ValueError(f'Reference walk index should be an integer >= 0 and < {G.num_walks}')

        G.add_reference_path(walks[reference_index])

        if nodeinfo_file:
            print("Reading nodeinfo file")
            G.read_nodeinfo(nodeinfo_file)
        else:
            print("Assigning node directions")
            assign_node_directions(G, G.reference_path)

        walk_start_nodes = [walk[0] for walk in walks]
        walk_end_nodes = [walk[-1] for walk in walks]
        G.add_terminal_nodes(walk_start_nodes=walk_start_nodes, walk_end_nodes=walk_end_nodes)

        if edgeinfo_file:
            print("Reading edgeinfo file")
            G.read_edgeinfo(edgeinfo_file)
        else:
            print("Computing reference tree")
            G.compute_edge_weights(walks)
            G.compute_reference_tree()
            print("Computing branch points")
            G.annotate_branch_points()

        if not nodeinfo_file:
            print("Computing positions")
            G.compute_binode_positions()

        if return_walks:
            return G, walks, walk_sample_names

        return G

    def __init__(self,
                 directed_graph: nx.classes.digraph.DiGraph = None,
                 reference_tree: nx.classes.digraph.DiGraph = None,
                 reference_path: list[str] = None,
                 variant_edges: set = None
                 ):
        if directed_graph is None:
            directed_graph = nx.DiGraph()
        if reference_tree is None:
            reference_tree = nx.DiGraph()
        super().__init__(directed_graph)
        self.reference_tree = reference_tree
        self.reference_path = reference_path if reference_path else []
        self.variant_edges = variant_edges if variant_edges else {}
        self.number_of_biedges = np.sum(
            [count_or_not for _, _, count_or_not in directed_graph.edges(data='is_representative')]
        )

    def variant_edges_summary(self):
        summary_dict = dict()
        for edge in tqdm(sorted(list(self.variant_edges))):
            if self.is_inversion(edge):
                summary_dict['inversion'] = summary_dict.get('inversion', 0) + 1
            if self.is_replacement(edge):
                summary_dict['replacement'] = summary_dict.get('replacement', 0) + 1
            if self.is_insertion(edge):
                summary_dict['insertion'] = summary_dict.get('insertion', 0) + 1
            if self.is_snp(edge):
                summary_dict['snps'] = summary_dict.get('snps', 0) + 1
            if self.is_mnp(edge):
                summary_dict['mnps'] = summary_dict.get('mnps', 0) + 1
            if self.is_crossing_edge(edge):
                summary_dict['crossing_edges'] = summary_dict.get('crossing_edges', 0) + 1
            if self.is_back_edge(edge):
                summary_dict['back_edges'] = summary_dict.get('back_edges', 0) + 1
            if self.is_forward_edge(edge):
                summary_dict['forward_edges'] = summary_dict.get('forward_edges', 0) + 1
        summary_dict['total'] = len(self.variant_edges)
        # Desired order of keys
        key_order = ['snps', 'mnps', 'insertion', 'replacement', 'inversion', 'crossing_edges', 'back_edges', 'forward_edges', 'total']
        # Creating a new dict with the desired order
        summary_dict = {key: summary_dict[key] for key in key_order}
        return summary_dict

    def is_inversion(self, edge):
        u, v = edge
        return self.nodes[u]['direction'] != self.nodes[v]['direction']

    def is_in_tree(self, edge):
        return self.representative_edge(edge) in self.reference_tree

    def is_back_edge(self, edge):
        if self.is_inversion(edge):
            return False
        if self.is_in_tree(edge):
            return False
        positive_edge = edge if self.nodes[edge[0]]['direction'] == 1 else _edge_complement(edge)
        branch_point = self.edges[positive_edge]['branch_point']
        return branch_point == positive_edge[1]

    def is_forward_edge(self, edge):
        if self.is_inversion(edge):
            return False
        if self.is_in_tree(edge):
            return False
        positive_edge = edge if self.nodes[edge[0]]['direction'] == 1 else _edge_complement(edge)
        branch_point = self.edges[positive_edge]['branch_point']
        return branch_point == positive_edge[0]

    def is_crossing_edge(self, edge):
        if self.is_inversion(edge):
            return False
        if self.is_in_tree(edge):
            return False
        if self.is_forward_edge(edge):
            return False
        if self.is_back_edge(edge):
            return False
        return True

    def is_insertion(self, edge):
        if not self.is_crossing_edge(edge):
            return False
        ref, _, _, _ = self.ref_alt_alleles(edge)
        return len(ref) == 0

    def is_replacement(self, edge):
        if not self.is_crossing_edge(edge):
            return False
        ref, alt, _, _ = self.ref_alt_alleles(edge)
        return len(ref) != 0 and len(alt) != 0 and len(ref) != len(alt)

    def is_snp(self, edge):
        if not self.is_crossing_edge(edge):
            return False
        ref, alt, _, _ = self.ref_alt_alleles(edge)
        return len(ref) == len(alt) == 1

    def is_mnp(self, edge):
        if not self.is_crossing_edge(edge):
            return False
        ref, alt, _, _ = self.ref_alt_alleles(edge)
        return len(ref) == len(alt) > 1

    def add_reference_path(self, reference_walk):
        self.reference_path = reference_walk
        for node in reference_walk:
            self.nodes[node]['on_reference_path'] = 1
            self.nodes[_node_complement(node)]['on_reference_path'] = 1


    def add_terminal_nodes(self, walk_start_nodes: list[str]=None, walk_end_nodes: list[str]=None):
        """Add two terminal binodes, +_terminus and -_terminus, to the graph. A valid walk proceeds
        from +_terminus_+ to -_terminus_+ or from -_terminus_- to +_terminus_-."""

        positive_subgraph = self.subgraph([n for n, direction in self.nodes(data="direction") if direction == 1])

        # source and sink nodes for positive-direction walks; complementary nodes are sink and source nodes for
        # negative-direction walks, respectively
        source_nodes = {node for node, degree in positive_subgraph.in_degree() if degree == 0}
        sink_nodes = {node for node, degree in positive_subgraph.out_degree() if degree == 0}

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

        plus_terminus, minus_terminus = self.termini
        self.add_binode(plus_terminus)
        self.add_binode(minus_terminus)

        for source_node in source_nodes:
            self.add_biedge(plus_terminus + '_+', source_node)

        for sink_node in sink_nodes:
            self.add_biedge(sink_node, minus_terminus + '_+')

        # reference path is assumed to be in the positive direction
        self.reference_path = [plus_terminus + '_+'] + self.reference_path + [minus_terminus + '_+']
        self.nodes[plus_terminus + '_+']['on_reference_path'] = 1
        self.nodes[plus_terminus + '_-']['on_reference_path'] = 1
        self.nodes[minus_terminus + '_+']['on_reference_path'] = 1
        self.nodes[minus_terminus + '_-']['on_reference_path'] = 1

        chromosome_length = self.nodes[self.reference_path[-2]]['position']
        self.nodes[plus_terminus + '_+']['position'] = 0
        self.nodes[plus_terminus + '_-']['position'] = 0
        self.nodes[minus_terminus + '_+']['position'] = chromosome_length
        self.nodes[minus_terminus + '_-']['position'] = chromosome_length
        self.nodes[plus_terminus + '_+']['forward_position'] = 0
        self.nodes[plus_terminus + '_-']['forward_position'] = 0
        self.nodes[minus_terminus + '_+']['forward_position'] = chromosome_length
        self.nodes[minus_terminus + '_-']['forward_position'] = chromosome_length

    def compute_reference_tree(self):

        # reference tree contains positive-direction nodes only, and no inversion edges
        positive_subgraph = self.subgraph([n for n, direction in self.nodes(data="direction") if direction == 1])
        self.reference_tree = max_weight_dfs_tree(positive_subgraph, reference_path=self.reference_path)

        for edge in positive_subgraph.edges():
            edge_in_tree = self.reference_tree.has_edge(*edge)
            self.edges[edge]['is_in_tree'] = edge_in_tree
            self.edges[_edge_complement(edge)]['is_in_tree'] = edge_in_tree

        # Variant edges are those not in the tree
        self.variant_edges = {(u, v) for u, v, data in self.edges(data=True)
                              if data['is_representative'] and not data['is_in_tree']}

    def write_vcf(self,
                  walks: list,
                  sample_names: list,
                  vcf_filename: str,
                  tree_filename: str,
                  chr_name: str,
                  size_threshold: int = 200,
                  walkup_limit: int = 50) -> None:
        # 'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'sample1', 'sample2', ...

        meta_info = f'##fileformat=VCFv4.2\n'
        meta_info += f'##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        meta_info += f'##INFO=<ID=NR,Number=1,Type=String,Description="The reference allele if it is not on reference path.">\n'
        meta_info += f'##INFO=<ID=PU,Number=1,Type=Integer,Description="Position of U (left node of variant edge)">\n'
        meta_info += f'##INFO=<ID=PV,Number=1,Type=Integer,Description="Position of V (right node of variant edge)">\n'
        meta_info += f'##INFO=<ID=LU,Number=1,Type=Integer,Description="Line number in the tree file">\n'
        meta_info += f'##INFO=<ID=LV,Number=1,Type=Integer,Description="Line number in the tree file">\n'
        meta_info += f'##INFO=<ID=RL,Number=1,Type=Integer,Description="Reference allele reach search limit">\n'
        meta_info += f'##INFO=<ID=AL,Number=1,Type=Integer,Description="Alternative allele reach search limit">\n'
        meta_info += f'##contig=<ID={chr_name[3:]}>\n'

        order: list = self.write_tree(tree_filename)
        node_to_line = {node: line for line, node in enumerate(order)}

        with open(vcf_filename, 'w') as file:
            # Write the header row
            sample_genotype_dict = self.integrate_genotype_by_sample(sample_names, walks)
            sample_genotype_dict.pop('GRCh38')
            sorted_sample_names = sorted(sample_genotype_dict.keys())
            header_names = list(self.vcf_attribute_names) + sorted_sample_names
            file.write(meta_info)
            file.write('#'+'\t'.join(header_names) + '\n')

            for u, v in tqdm(sorted(list(self.variant_edges), key=lambda x: self.nodes[x[0]]['position'])):
                if self.nodes[u]['direction'] != self.nodes[v]['direction']:
                    continue  # TODO how to handle inversions?

                if self.nodes[u]['direction'] == -1:
                    u, v = _edge_complement((u,v))
                edge = (u, v)

                ref_allele, alt_allele, last_letter_of_branch_point, branch_point, ref_limit, alt_limit = self.ref_alt_alleles(edge,
                                                                                                         walkup_limit=walkup_limit,
                                                                                                         return_search_bool=True)
                # ref_allele, alt_allele, last_letter_of_branch_point, branch_point = '', '', '', ''

                if len(ref_allele) == 0 or len(alt_allele) == 0:
                    ref_allele = last_letter_of_branch_point + ref_allele
                    alt_allele = last_letter_of_branch_point + alt_allele

                ref = ref_allele[:size_threshold]
                alt = alt_allele[:size_threshold]

                allele_data_list = []

                if self.nodes[v]['on_reference_path'] == 1:
                    new_ref = '.'
                else:
                    new_ref = ref
                    ref = '.'

                # 'CHROM'
                allele_data_list.append(chr_name)
                # 'POS'
                if (self.is_snp(edge) or self.is_mnp(edge)) and self.nodes[v]['on_reference_path'] == 1:
                    allele_data_list.append(str(int(self.nodes[u]['position']) + 1))
                else:
                    allele_data_list.append(str(int(self.nodes[u]['position'])))
                # 'ID'
                allele_data_list.append('.')
                # 'REF'
                allele_data_list.append(ref)
                # 'ALT'
                allele_data_list.append(alt)
                # 'QUAL'
                allele_data_list.append('60')
                # 'FILTER'
                allele_data_list.append('PASS')
                # 'INFO'
                if u not in node_to_line:
                    print(u, list(node_to_line.items()))

                allele_data_list.append(f'NR={new_ref};'
                                        f'PU={int(self.nodes[u]["position"])};'
                                        f'PV={int(self.nodes[v]["position"])};'
                                        f'LU={int(node_to_line[u])};'
                                        f'LV={int(node_to_line[v])};'
                                        f'RL={int(ref_limit)};'
                                        f'AL={int(alt_limit)};')
                # 'FORMAT'
                allele_data_list.append('GT')

                for sample_name in sorted_sample_names:
                    if edge in sample_genotype_dict[sample_name]:
                        if str(sample_name).startswith("CHM"):
                            allele_data_list.append(str(sample_genotype_dict[sample_name][edge][:1]))
                        else:
                            allele_data_list.append(str(sample_genotype_dict[sample_name][edge]))
                    else:
                        if str(sample_name).startswith("CHM"):
                            allele_data_list.append('0')
                        else:
                            allele_data_list.append('0|0')

                file.write('\t'.join(allele_data_list) + '\n')

    def write_variant_node_ids(self, filename: str) -> None:
        with open(filename, 'w') as file:
            file.write("Node_u_id,Node_v_id\n")
            for u, v in self.variant_edges:
                file.write(f"{u},{v}\n")

    def write_tree(self, filename: str) -> list:
        with open(filename, 'w') as file:
            line_in_file = {}
            current_position = 0

            def write_node(node):
                if self.reference_tree.in_degree[node] == 0:
                    parent_position_in_file = -1
                else:
                    parent = next(self.reference_tree.predecessors(node))
                    parent_position_in_file = line_in_file[parent] if parent else -1
                node_sequence = self.nodes[node]['sequence']
                file.write(f"{node},{parent_position_in_file},{node_sequence}\n")

            order = list(nx.topological_sort(self.reference_tree))
            for line, u in enumerate(order):
                write_node(u)
                line_in_file[u] = line

            return order

    def write_edgeinfo(self, filename: str) -> None:
        def to_string(edge_attribute):
            if type(edge_attribute) is bool:
                return '1' if edge_attribute else '0'
            return str(edge_attribute)

        with open(filename, 'w') as file:
            # Write the header row
            file.write(','.join(self.biedge_attribute_names) + '\n')

            # Sort the edges by 'index' attribute, restricting to representatives
            sorted_edges = self.sorted_biedge_representatives

            for n, edge in enumerate(sorted_edges):
                _, _, data = edge
                assert n == data['index'], 'Something is wrong with edge indices'
                # Write the edge information to the file; map True->'1', False->'0'
                edge_data_list = [to_string(data[key]) for key in self.biedge_attribute_names]
                file.write(','.join(edge_data_list) + '\n')

    def write_nodeinfo(self, filename: str) -> None:

        attributes = [attribute for attribute in self.node_attribute_names if attribute != 'sequence']
        with open(filename, 'w') as file:
            file.write('node,')
            file.write(','.join(attributes) + '\n')

            for node, data in self.nodes(data=True):
                if self.is_terminal_node(node):
                    continue
                file.write(f'{node},')
                file.write(','.join([str(data[key]) for key in attributes]) + '\n')

    def read_edgeinfo(self, filename: str) -> None:

        def from_string(s: str):
            try:
                float(s)
                return float(s)
            except ValueError:
                return s

        sorted_edges = self.sorted_biedge_representatives

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
                    self.edges[edge][key] = from_string(parts[i])
                    self.edges[_edge_complement(edge)][key] = from_string(parts[i])

        # Define variant edge representatives
        self.variant_edges = {(u, v) for u, v, data in self.edges(data=True)
                              if data['is_representative'] and not data['is_in_tree']}

        # Define reference tree
        for u, v, is_in_tree in self.edges(data='is_in_tree'):
            if not is_in_tree:
                continue

            # Only include forward direction nodes
            if self.nodes[u]['direction'] == 1:
                self.reference_tree.add_edge(u, v)

        # Connect non-terminus roots of the reference tree (forest) with terminus
        root = self.termini[0] + '_+'
        source_nodes = [node for node, in_degree in self.reference_tree.in_degree if in_degree == 0]
        for source_node in source_nodes:
            if source_node == root:
                continue
            assert not self.reference_tree.has_edge(root, source_node)
            self.reference_tree.add_edge(root, source_node)

    def read_nodeinfo(self, filename: str) -> None:
        def from_string(s: str):
            try:
                float(s)
                return float(s)
            except ValueError:
                return s

        with open(filename, 'r') as file:
            line = next(file).strip().split(',')
            attributes = {key: idx for idx, key in enumerate(line) if key != 'node'}
            for line in file:
                parts = line.strip().split(',')
                node = parts[0]
                for key, idx in attributes.items():
                    self.nodes[node][key] = from_string(parts[idx])



    def add_binode(self, binode: str, seq: str = ''):
        node_data = {key: 0 for key in self.node_attribute_names}
        node_data['sequence'] = seq
        node_data['direction'] = 1

        self.add_node(str(binode) + '_+', **node_data)

        node_data['sequence'] = _sequence_reversed_complement(seq)
        node_data['direction'] = -1

        self.add_node(str(binode) + '_-', **node_data)

    def add_biedge(self, node1: str, node2: str, *, weight: int = 0):

        if self.has_edge(node1, node2):
            raise ValueError(f'Attempted to add duplicate biedge: {node1}, {node2}')

        edge_data = {key: 0 for key in self.biedge_attribute_names}
        edge_data['weight'] = weight
        edge_data['index'] = self.number_of_biedges
        edge_data['is_representative'] = True
        self.add_edge(node1, node2, **edge_data)

        edge_data['is_representative'] = False
        self.add_edge(_node_complement(node2), _node_complement(node1), **edge_data)
        self.number_of_biedges += 1

    def representative_edge(self, edge: tuple):
        return edge if self.edges[edge]['is_representative'] else _edge_complement(edge)

    def positive_variant_edge(self, edge: tuple):
        return edge if self.nodes[edge[0]]['direction'] == 1 or self.is_inversion(edge) else _edge_complement(edge)
    def compute_edge_weights(self, walks):
        for walk in walks:

            for u, v in zip(walk[:-1], walk[1:]):
                self.edges[u, v]['weight'] += 1
                self.edges[_node_complement(v), _node_complement(u)]['weight'] += 1

    def annotate_branch_points(self) -> None:
        """
        Computes the branch point, i.e. the lowest common ancestor in the reference tree, of each variant edge.
        """

        def positive_node(node):
            return node if self.nodes[node]['direction'] == 1 else _node_complement(node)

        positive_direction_variants = [(positive_node(u), positive_node(v)) for u,v in self.variant_edges]

        # For each variant edge (u,v), yield the lowest common ancestor of u and v in the tree
        branch_point_tuples = nx.tree_all_pairs_lowest_common_ancestor(
            self.reference_tree,
            root=self.termini[0]+'_+',
            pairs=positive_direction_variants
        )
        branch_points = dict(branch_point_tuples)

        for edge in self.variant_edges:
            branch_point = branch_points[positive_node(edge[0]), positive_node(edge[1])]
            self.edges[edge]['branch_point'] = branch_point
            self.edges[_edge_complement(edge)]['branch_point'] = branch_point

    def walk_up_tree(self, ancestor, descendant, search_limit=50) -> tuple[list, bool]:
        u = descendant
        result = []
        search_count = 0
        reach_limit = False
        while u != ancestor and search_count < search_limit:
            search_count += 1
            result.append(u)
            u = next(self.reference_tree.predecessors(u))
        result.append(ancestor)
        if (search_count + 1) == search_limit:
            reach_limit = True
        return result, reach_limit

    # TODO think about strandedness and desired behavior; currently this maps [-, -] variant edges to [+, +] silently
    def ref_alt_alleles(self, variant_edge: tuple,
                        walkup_limit: int = 50,
                        return_search_bool: bool=False,
                        ):
        """
        Computes the reference allele and alternative allele of the branch point for each variant edge.
        :param variant_edges: list of tuples (u,v).
        :return: dict mapping variant edges to tuples (ref, alt)
        """

        variant_edge = self.representative_edge(variant_edge)
        if self.is_in_tree(variant_edge):
            raise ValueError("Ref and alt alleles are only defined for variant edges")

        u, v = variant_edge
        branch_point = self.edges[u, v]['branch_point']
        if self.nodes[u]['direction'] == -1:
             u, v = _edge_complement((u, v))

        # inversion edge
        if self.nodes[v]['direction'] == -1:
            v = _node_complement(v)

        ref_path, ref_search_limit = self.walk_up_tree(branch_point, v, walkup_limit)
        alt_path, alt_search_limit = self.walk_up_tree(branch_point, u, walkup_limit)

        # alt allele sometimes includes and sometime excludes the branch point and u, depending on edge type
        is_back_edge = branch_point == v
        is_forward_edge = branch_point == u
        if is_back_edge:
            alt_path = alt_path[::-1]
        elif is_forward_edge:
            alt_path = []
        else:
            alt_path = alt_path[-2::-1]

        # ref path always excludes both the branch point and v
        ref_path = ref_path[-2:0:-1]

        alt_allele = ''
        if alt_search_limit:
            alt_allele = '.'
        else:
            for node in alt_path:
                alt_allele += self.nodes[node]['sequence']

        ref_allele = ''
        if ref_search_limit:
            ref_allele = '.'
        else:
            for node in ref_path:
                ref_allele += self.nodes[node]['sequence']

        branch_sequence = self.nodes[branch_point]['sequence']
        if not branch_sequence:
            branch_sequence = 'N'
        last_letter_of_branch_point = branch_sequence[-1]

        if return_search_bool:
            return ref_allele, alt_allele, last_letter_of_branch_point, branch_point, ref_search_limit, alt_search_limit
        else:
            return ref_allele, alt_allele, last_letter_of_branch_point, branch_point

    def integrate_genotype_by_sample(self, sample_names, walks):
        sample_haplotype_dict = defaultdict(Counter)
        sample_diploid_dict = defaultdict(list)

        for sample_name, walk in zip(sample_names, walks):
            walk_haplotype_dict = self.genotype(walk)
            sample_haplotype_dict[sample_name].update(walk_haplotype_dict)

        for k, v in sample_haplotype_dict.items():
            sample_name, _ = k.split('_')
            sample_diploid_dict[sample_name].append(v)

        def concatenate_dicts(dicts: list[Counter]):
            result_dict = dict()

            dict_1 = dicts[0]

            if len(dicts) == 1:
                dict_2 = dict()
            elif len(dicts) == 2:
                dict_2 = dicts[1]
            else:
                raise ValueError("This should never happen")

            all_keys = set(dict_1.keys()).union(set(dict_2.keys()))

            for key in all_keys:
                if key in dict_1 and key in dict_2:
                    result_dict[key] = f"{dict_1[key]}|{dict_2[key]}"
                elif key in dict_1:
                    result_dict[key] = f"{dict_1[key]}|0"
                elif key in dict_2:
                    result_dict[key] = f"0|{dict_2[key]}"
                else:
                    raise ValueError("This should never happen")

            return result_dict

        return {k: concatenate_dicts(v) for k, v in sample_diploid_dict.items()}


    def genotype(self, walk: list[str]) -> dict:

        # Append start and end nodes to walk
        start = [self.termini[0] + '_+' if self.nodes[walk[0]]['direction'] == 1 else self.termini[1] + '_-']
        end = [self.termini[1] + '_+' if self.nodes[walk[-1]]['direction'] == 1 else self.termini[0] + '_-']
        walk = start + walk + end

        genotype = {}
        for e in zip(walk[:-1], walk[1:]):
            if not self.has_edge(*e):
                raise ValueError(f"Specified list contains edge {e} which is not present in the graph")

            if not self.edges[e]['is_representative']:
                e = _edge_complement(e)

            if self.edges[e]['is_in_tree']:
                continue

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

    def compute_binode_positions(self):
        for node in self.reference_tree.nodes():
            self.nodes[node]['position'] = -inf
            self.nodes[node]['forward_position'] = inf

        current_position = 0
        for u in self.reference_path:
            current_position += len(self.nodes[u]['sequence'])
            self.nodes[u]['position'] = current_position
            self.nodes[u]['forward_position'] = current_position

        order = list(nx.topological_sort(self.reference_tree))
        for u in order[1:]: # skip the root
            predecessor = next(self.reference_tree.predecessors(u))
            self.nodes[u]['position'] = np.maximum(self.nodes[u]['position'],
                                               self.nodes[predecessor]['position'])
            self.nodes[_node_complement(u)]['position'] = self.nodes[u]['position']

        for u in reversed(order[:-1]):
            successors = [v for _, v, is_back_edge in self.out_edges(u, data='is_back_edge') if not is_back_edge]
            if successors:
                successors_minimum_position = np.min([self.nodes[v]['forward_position'] for v in successors])
            else:
                successors_minimum_position = inf

            self.nodes[u]['forward_position'] = np.minimum(self.nodes[u]['forward_position'],
                                                            successors_minimum_position)
            self.nodes[_node_complement(u)]['forward_position'] = self.nodes[u]['forward_position']

    def get_variants_at_interval(self, half_open_interval: tuple[int, int], exclude_root_edges=True) -> list:
        start, end = half_open_interval
        result = []

        for edge in self.variant_edges:
            u, v = edge
            positions = [self.nodes[u]['position'], self.nodes[v]['position']]
            x, y = min(positions), max(positions)
            if x < end and y >= start:
                result.append(edge)

        if exclude_root_edges:
            result = [(u,v) for u,v in result if u != '+_terminus_+' and v != '+_terminus_-']

        return result

    def classify_triallelic_bubble(self, endpoint_binodes: list, variants: list) -> str:
        """
        Classifies a top-level triallelic superbubble containing 2 variants as either:
        overlapping: one pair of the three alleles shares one subsequence
        interlocking: two different pairs of alleles share one subsequence each
        nested: one pair of alleles shares two noncontiguous subsequences
        properly_triallelic: no pair of alleles shares a subsequence
        # TODO unsure if this handles inversions properly

        :param endpoint_binodes: starting and ending binodes of the superbubble
        :param variants: list of variants within the superbubble; there must be two of them
        :return: the class of superbubble
        """
        if len(variants) != 2:
            raise ValueError

        endpoint_nodes = [node + '_+' for node in endpoint_binodes]
        endpoint_nodes += [node + '_-' for node in endpoint_binodes]
        if not all([self.has_node(node) for node in endpoint_nodes]):
            raise ValueError

        # Detect which node of each binode demarcates the superbubble
        variant_branch_points = [self.edges[e]['branch_point'] for e in variants]
        start_nodes = [node for node in endpoint_nodes if node in variant_branch_points]
        variant_end_points = [v for _, v in variants]
        end_nodes = [node for node in endpoint_nodes if node in variant_end_points]
        assert len(start_nodes) == 1 and len(end_nodes) == 1, \
            f"Found {len(start_nodes)} possible start nodes and {len(end_nodes)} possible end nodes"

        start_node = start_nodes[0]
        end_node = end_nodes[0]
        start_degree = self.out_degree(start_node)
        end_degree = self.in_degree(end_node)
        assert start_degree <= 3 and end_degree <= 3, \
            f"Starting and ending nodes ({start_node} and {end_node}) of the bubble had degree {start_degree} and {end_degree}"

        if start_degree == 3 and end_degree == 3:
            return 'properly_triallelic'

        if start_degree == 3 or end_degree == 3:
            return 'overlapping'

        # e.g., [(0,1), (0,2), (1,2), (1,3), (2,3)] with variant edges [(1,2), (1,3)]
        if all([node == start_node for node in variant_branch_points]):
            return 'interlocking'

        for branch_point, end_point in zip(variant_branch_points, variant_end_points):
            if branch_point == start_node and end_point == end_node:
                return 'nested'

        # e.g., [(0,1), (0,2), (1,2), (1,3), (2,3)] with variant edges [(0,2), (1,3)]
        return 'interlocking'





