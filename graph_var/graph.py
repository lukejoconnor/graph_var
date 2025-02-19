from functools import lru_cache
from linecache import cache
from math import inf
import networkx as nx
import numpy as np
from .utils import read_gfa, node_complement, edge_complement, sequence_complement, walk_complement, group_walks_by_name, nearly_identical_alleles
from .search_tree import assign_node_directions, max_weight_dfs_tree
import os
from collections import defaultdict, Counter
from tqdm import tqdm
from typing import Any, Union, Optional

class PangenomeGraph(nx.DiGraph):
    reference_tree: nx.classes.digraph.DiGraph
    reference_path: list[str]
    variant_edges: set
    number_of_biedges: int  # Needed because nx.number_of_edges() runs in O(number of edges)

    # A valid walk proceeds from either +_terminus_+ to -_terminus_+ or from -_terminus_- to +_terminus_-
    @property
    def termini(self) -> tuple[str, str]:
        return '+_terminus', '-_terminus'

    def is_terminal(self, node_or_edge) -> bool:
        if isinstance(node_or_edge, str):
            return node_or_edge[:-2] in self.termini
        elif isinstance(node_or_edge, tuple):
            return self.is_terminal(node_or_edge[0]) or self.is_terminal(node_or_edge[1])
        else:
            raise TypeError

    # Each biedge has a representative edge, whichever is in the .gfa file
    @property
    def sorted_biedge_representatives(self) -> list[str]:
        edges = [edge_with_data for edge_with_data in self.edges(data=True) if edge_with_data[2]['is_representative']]
        return sorted(edges, key=lambda edge: edge[2]['index'])

    def sorted_variant_edges(self, exclude_terminus=True) -> list[str]:
        
        if exclude_terminus:
            sorted_vars = [edge for edge in self.variant_edges if not self.is_terminal(edge)]
        else:
            sorted_vars = self.variant_edges

        sorted_vars = sorted(sorted_vars, key=lambda x:
                      (self.nodes[self.positive_variant_edge(x)[0]]['position'],
                       int(self.nodes[self.positive_variant_edge(x)[0]]["distance_from_reference"])))
        return sorted_vars

    @property
    def biedge_attribute_names(self) -> tuple:
        return 'index', 'weight', 'is_in_tree', 'branch_point', 'is_back_edge'

    @property
    def node_attribute_names(self) -> tuple:
        return 'direction', 'sequence', 'position', 'right_position', 'distance_from_reference', 'on_reference_path'

    @property
    def vcf_attribute_names(self) -> tuple:
        return 'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'

    @property
    def number_of_binodes(self) -> int:
        return self.number_of_nodes() // 2

    def on_reference_path(self, node_or_edge):
        if type(node_or_edge) is tuple:
            if not self.has_edge(node_or_edge):
                raise ValueError("Graph does not have edge {node_or_edge}")
            return self.nodes[node_or_edge[1]]['on_reference_path']
        elif type(node_or_edge) is str:
            if not self.has_node(node_or_edge):
                raise ValueError("Graph does not have node {node_or_edge}")
            return self.nodes[node_or_edge]['on_reference_path']

    def parent_in_tree(self, node: str) -> str:
        if self.reference_tree.in_degree(node) == 0:
            return None
        return next(self.reference_tree.predecessors(node))

    def position(self, node_or_edge: Any) -> int:
        if isinstance(node_or_edge, str):
            return self.nodes[node_or_edge]['position']
        elif isinstance(node_or_edge, tuple):
            return self.position(node_or_edge[0]), self.position(node_or_edge[1])
        else:
            raise TypeError

    def right_position(self, node_or_edge: Any) -> int:
        if isinstance(node_or_edge, str):
            return self.nodes[node_or_edge]['right_position']
        elif isinstance(node_or_edge, tuple):
            return self.right_position(node_or_edge[0]), self.right_position(node_or_edge[1])
        else:
            raise TypeError

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
        :param gfa_file: path to a file name ending in .gfa
        :param reference_path_index: one of the walks in the gfa file can be specified as the reference path; this
        is the index of that walk. If unspecified, the first walk whose name is 'GRCh38' is used, if any.
        :param edgeinfo_file: path to a previously-computed .edgeinfo file, to avoid re-doing work
        :param nodeinfo_file: path to a previously-computed .nodeinfo file, to avoid re-doing work
        :param return_walks: if True, return a tuple (graph_object, walks, walk_names) where the walks are
        read from the .gfa file; if False, return graph_object
        :param compressed: set to True in order to read a gzipped .gfa file
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
            raise ValueError(f'Reference walk index should be an integer >= 0')

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
            G.compute_binode_right_positions()

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

    def identify_variant_type(self, edge: tuple[str, str]) -> str:
        var_type = None
        if self.is_inversion(edge):
            if var_type is not None:
                raise KeyError(f'Variant, {edge}, has dual type.')
            var_type = 'inversion'
        if self.is_replacement(edge):
            if var_type is not None:
                raise KeyError(f'Variant, {edge}, has dual type.')
            var_type = 'replacement'
        if self.is_insertion(edge):
            if var_type is not None:
                raise KeyError(f'Variant, {edge}, has dual type.')
            var_type = 'insertion'
        if self.is_snp(edge):
            if var_type is not None:
                raise KeyError(f'Variant, {edge}, has dual type.')
            var_type = 'snp'
        if self.is_mnp(edge):
            if var_type is not None:
                raise KeyError(f'Variant, {edge}, has dual type.')
            var_type = 'mnp'
        if self.is_back_edge(edge):
            if var_type is not None:
                raise KeyError(f'Variant, {edge}, has dual type.')
            var_type = 'repeat'
        if self.is_forward_edge(edge):
            if var_type is not None:
                raise KeyError(f'Variant, {edge}, has dual type.')
            var_type = 'deletion'
        return var_type

    def is_inversion(self, edge: tuple[str, str]) -> bool:
        u, v = edge
        return self.direction(u) != self.direction(v)

    def is_in_tree(self, edge: tuple[str, str]) -> bool:
        return self.representative_edge(edge) in self.reference_tree

    def is_back_edge(self, edge: tuple[str, str]) -> bool:
        if self.is_inversion(edge):
            return False
        return self.edges[edge]['is_back_edge']

    def is_forward_edge(self, edge: tuple[str, str]) -> bool:
        if self.is_inversion(edge):
            return False
        if self.is_in_tree(edge):
            return False
        positive_edge = edge if self.direction(edge[0]) == 1 else edge_complement(edge)
        branch_point = self.edges[positive_edge]['branch_point']
        return branch_point == positive_edge[0]

    def is_crossing_edge(self, edge: tuple[str, str]) -> bool:
        if self.is_inversion(edge):
            return False
        if self.is_in_tree(edge):
            return False
        if self.is_forward_edge(edge):
            return False
        if self.is_back_edge(edge):
            return False
        return True

    def is_insertion(self, edge: tuple[str, str]) -> bool:
        if not self.is_crossing_edge(edge):
            return False
        ref, _, _, _ = self.ref_alt_alleles(edge)
        return len(ref) == 0

    def is_replacement(self, edge: tuple[str, str]) -> bool:
        if not self.is_crossing_edge(edge):
            return False
        ref, alt, _, _ = self.ref_alt_alleles(edge)
        return len(ref) != 0 and len(alt) != 0 and len(ref) != len(alt)

    def is_snp(self, edge: tuple[str, str]) -> bool:
        if not self.is_crossing_edge(edge):
            return False
        ref, alt, _, _ = self.ref_alt_alleles(edge)
        return len(ref) == len(alt) == 1

    def is_mnp(self, edge: tuple[str, str]) -> bool:
        if not self.is_crossing_edge(edge):
            return False
        ref, alt, _, _ = self.ref_alt_alleles(edge)
        return len(ref) == len(alt) > 1

    def add_reference_path(self, reference_walk):
        self.reference_path = reference_walk
        for node in reference_walk:
            self.nodes[node]['on_reference_path'] = 1
            self.nodes[node_complement(node)]['on_reference_path'] = 1


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
                {node for node in walk_start_nodes if self.direction(node) == 1}
            )
            sink_nodes = sink_nodes.union(
                {node_complement(node) for node in walk_start_nodes if self.direction(node) == -1}
            )
        if walk_end_nodes:
            sink_nodes = sink_nodes.union(
                {node for node in walk_end_nodes if self.direction(node) == 1}
            )
            source_nodes = source_nodes.union(
                {node_complement(node) for node in walk_end_nodes if self.direction(node) == -1}
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

    def compute_reference_tree(self):
        """
        Computes the reference tree, a DFS spanning tree of the positively-oriented subgraph; defines variant edges
        as those that are not in the reference tree or its complement.
        """
        # reference tree contains positive-direction nodes only, and no inversion edges
        positive_subgraph = self.subgraph([n for n, direction in self.nodes(data="direction") if direction == 1])
        self.reference_tree = max_weight_dfs_tree(positive_subgraph, reference_path=self.reference_path)

        for edge in positive_subgraph.edges():
            edge_in_tree = self.reference_tree.has_edge(*edge)
            self.edges[edge]['is_in_tree'] = edge_in_tree
            self.edges[edge_complement(edge)]['is_in_tree'] = edge_in_tree

        # Variant edges are those not in the tree
        self.variant_edges = {(u, v) for u, v, data in self.edges(data=True)
                              if data['is_representative'] and not data['is_in_tree']}

    def write_vcf(self,
                  walks: list,
                  sample_names: list,
                  vcf_filename: str,
                  tree_filename: str,
                  chr_name: str,
                  size_threshold: int = None,
                  walkup_limit: int = inf,
                  exclude_terminus: bool = True,
                  check_degenerate: bool = False,
                  ) -> None:
        """
        Writes the variant call format (vcf) file and the tree file. The vcf file contains the reference and
        alternative alleles, and the tree file contains the node names and corresponding sequences.
        :param walks: the list of walks extracted from the .gfa file
        :param sample_names: the haplotype id from the .gfa file
        :param vcf_filename: the output vcf file path
        :param tree_filename: the output tree file path which contains the node and corresponding sequence
        :param chr_name: the chromosome name in the first column of output vcf file
        :param size_threshold: the truncation length of ref and alt sequence
        :param walkup_limit: the maximum length of the walk up to the branch point
        :param exclude_terminus: whether to exclude the terminus nodes
        :param check_degenerate: whether to exclude variants whose ref and alt alleles are identical
        :return:
        """
        # 'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'sample1', 'sample2', ...
        meta_info = f'##fileformat=VCFv4.2\n'
        meta_info += f'##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        meta_info += f'##INFO=<ID=OR,Number=1,Type=String,Description="Off-linear reference.">\n'
        meta_info += f'##INFO=<ID=VT,Number=1,Type=Integer,Description="Variant type.">\n'
        meta_info += f'##INFO=<ID=DR,Number=1,Type=Integer,Description="Distance of node u, v from reference.">\n'
        meta_info += f'##INFO=<ID=RC,Number=1,Type=String,Description="The REF allele count.">\n'
        meta_info += f'##INFO=<ID=AC,Number=1,Type=String,Description="The ALT allele count.">\n'
        meta_info += f'##INFO=<ID=CR,Number=1,Type=String,Description="The REF count for haplotype.">\n'
        meta_info += f'##INFO=<ID=CA,Number=1,Type=String,Description="The ALT count for haplotype.">\n'
        meta_info += f'##INFO=<ID=AN,Number=1,Type=String,Description="The non-missing value count.">\n'
        meta_info += f'##INFO=<ID=PU,Number=1,Type=Integer,Description="Position of U (left node of variant edge)">\n'
        meta_info += f'##INFO=<ID=PV,Number=1,Type=Integer,Description="Position of V (right node of variant edge)">\n'
        meta_info += f'##INFO=<ID=LU,Number=1,Type=Integer,Description="Line number in the tree file">\n'
        meta_info += f'##INFO=<ID=LV,Number=1,Type=Integer,Description="Line number in the tree file">\n'
        meta_info += f'##INFO=<ID=RL,Number=1,Type=Integer,Description="Reference allele reach search limit">\n'
        meta_info += f'##INFO=<ID=AL,Number=1,Type=Integer,Description="Alternative allele reach search limit">\n'
        meta_info += '##INFO=<ID=NIA,Number=0,Type=Flag,Description="Nearly identical alleles">\n'
        meta_info += f'##contig=<ID={chr_name[3:]}>\n'

        order: list = self.write_tree(tree_filename)
        node_to_line = {node: line for line, node in enumerate(order)}
        allele_count_dict = self.allele_count()

        with open(vcf_filename, 'w') as file:
            # Write the header row
            sample_walks_dict = group_walks_by_name(walks, sample_names)
            sample_data_dict = {sample_name: self.genotype_and_linear_coverage_by_sample(walks) for sample_name, walks in sample_walks_dict.items()}
            sample_missing_dict = {sample_name: set(self.get_missing_variants(data[2])) for sample_name, data in sample_data_dict.items()}

            #subject_ids = sorted({sample_name.split('_')[0] for sample_name in sample_walks_dict.keys() if not sample_name.startswith("GRCh38")})
            subject_ids = sorted({sample_name.split('_')[0] for sample_name in sample_walks_dict.keys()})

            header_names = list(self.vcf_attribute_names) + subject_ids
            file.write(meta_info)
            file.write('#'+'\t'.join(header_names) + '\n')

            for u, v in tqdm(self.sorted_variant_edges(exclude_terminus=exclude_terminus)):
                representative_variant_edge = (u, v)
                representative_ref_edge = self.representative_edge(self.reference_tree_edge(representative_variant_edge))

                if self.direction(u) == -1 and self.direction(v) == -1:
                    u, v = edge_complement((u, v))
                edge = (u, v)

                ref_allele, alt_allele, last_letter_of_branch_point, branch_point, ref_limit, alt_limit = self.ref_alt_alleles(edge,
                                                                                                         walkup_limit=walkup_limit,
                                                                                                         return_search_bool=True)

                if len(ref_allele) == 0 or len(alt_allele) == 0:
                    ref_allele = last_letter_of_branch_point + ref_allele
                    alt_allele = last_letter_of_branch_point + alt_allele

                if size_threshold:
                    ref = ref_allele[:size_threshold]
                    alt = alt_allele[:size_threshold]
                else:
                    ref = ref_allele
                    alt = alt_allele

                if check_degenerate:
                    if ref_allele == alt_allele:
                        continue

                allele_data_list = []

                if self.nodes[v]['on_reference_path'] == 1:
                    new_ref = '.'
                else:
                    new_ref = ref
                    ref = '.'

                # 'CHROM' 0
                allele_data_list.append(chr_name)
                # 'POS' 1
                allele_data_list.append(str(self.get_variant_position(edge)))
                # 'ID' 2
                allele_data_list.append(f"{representative_variant_edge}")
                # 'REF' 3
                allele_data_list.append(ref)
                # 'ALT' 4
                allele_data_list.append(alt)
                # 'QUAL' 5
                allele_data_list.append('60')
                # 'FILTER' 6
                allele_data_list.append('PASS')
                # 'INFO' 7
                allele_data_list.append(None)
                # 'FORMAT' 8
                allele_data_list.append('GT:CR:CA')

                # 'sample1', 'sample2', ... 9 - end
                AN = 0
                for sample_name in subject_ids:
                    if sample_name.startswith(("CHM", "GRCh")):
                        haplotype_name = sample_name+'_0'
                        counts = (f"{int(bool(sample_data_dict[haplotype_name][1].get(representative_variant_edge, 0)))}:"
                                  f"{sample_data_dict[haplotype_name][0].get(representative_ref_edge, 0)}:"
                                  f"{sample_data_dict[haplotype_name][1].get(representative_variant_edge, 0)}"
                        if representative_variant_edge not in sample_missing_dict[haplotype_name] else '.:.:.')
                        allele_data_list.append(counts)
                    else:
                        haplotype1_name = sample_name + '_1'
                        haplotype2_name = sample_name + '_2'
                        count_1 = (f"{int(bool(sample_data_dict[haplotype1_name][1].get(representative_variant_edge, 0)))}:"
                                  f"{sample_data_dict[haplotype1_name][0].get(representative_ref_edge, 0)}:"
                                  f"{sample_data_dict[haplotype1_name][1].get(representative_variant_edge, 0)}"
                                  if representative_variant_edge not in sample_missing_dict[haplotype1_name] else '.:.:.')
                        count_2 = (f"{int(bool(sample_data_dict[haplotype2_name][1].get(representative_variant_edge, 0)))}:"
                                   f"{sample_data_dict[haplotype2_name][0].get(representative_ref_edge, 0)}:"
                                   f"{sample_data_dict[haplotype2_name][1].get(representative_variant_edge, 0)}"
                                   if representative_variant_edge not in sample_missing_dict[haplotype2_name] else '.:.:.')
                        counts = f"{count_1}|{count_2}"
                        allele_data_list.append(counts)

                    for count in counts.split('|'):
                        if count != '.:.:.':
                            AN += 1

                INFO = (f'OR={new_ref};'
                        f'VT={self.identify_variant_type(edge)};'
                        f'DR={int(self.nodes[u]["distance_from_reference"])},{int(self.nodes[v]["distance_from_reference"])};'
                        f'RC={allele_count_dict[representative_variant_edge][0]};'
                        f'AC={allele_count_dict[representative_variant_edge][1]};'
                        f'AN={AN};'
                        f'PU={int(self.nodes[u]["position"])};'
                        f'PV={int(self.nodes[v]["position"])};'
                        f'LU={int(node_to_line[self.positive_node(u)])};'
                        f'LV={int(node_to_line[self.positive_node(v)])};'
                        f'RL={int(ref_limit)};'
                        f'AL={int(alt_limit)}')

                if nearly_identical_alleles(ref_allele, alt_allele):
                    INFO += ';NIA=1'

                allele_data_list[7] = INFO

                file.write('\t'.join(allele_data_list) + '\n')

    def write_tree(self, filename: str) -> list:
        """
        Writes the reference tree in a format such that it can be traversed without being loaded into memory.
        :param filename: save file path, ending by convention in .tree
        :return: order in which nodes were recorded in the file
        """
        with open(filename, 'w') as file:
            line_in_file = {}

            def write_node(node):
                if self.reference_tree.in_degree[node] == 0:
                    parent_position_in_file = -1
                else:
                    parent = self.parent_in_tree(node)
                    parent_position_in_file = line_in_file[parent] if parent else -1
                node_sequence = self.nodes[node]['sequence']
                file.write(f"{node},{parent_position_in_file},{node_sequence}\n")

            order = list(nx.topological_sort(self.reference_tree))
            for line, u in enumerate(order):
                write_node(u)
                line_in_file[u] = line

            return order

    def write_edgeinfo(self, filename: str) -> None:
        """
        Writes information computed about each edge such that it can be retrieved without re-computing it.
        :param filename: save file path, ending by convention in .edgeinfo
        """
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
        """
        Writes information about each node such that it can be retrieved without re-computing it.
        :param filename: save file path, ending by convention in .nodeinfo
        """
        attributes = [attribute for attribute in self.node_attribute_names if attribute != 'sequence']
        with open(filename, 'w') as file:
            file.write('node,')
            file.write(','.join(attributes) + '\n')

            for node, data in self.nodes(data=True):
                if self.is_terminal(node):
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
                    self.edges[edge_complement(edge)][key] = from_string(parts[i])

        # Define variant edge representatives
        self.variant_edges = {(u, v) for u, v, data in self.edges(data=True)
                              if data['is_representative'] and not data['is_in_tree']}

        # Define reference tree
        for u, v, is_in_tree in self.edges(data='is_in_tree'):
            if not is_in_tree:
                continue

            # Only include forward direction nodes
            if self.direction(u) == 1:
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
        """
        Adds a binode to the bidirected graph, comprising a node and its complement.
        """
        node_data = {key: 0 for key in self.node_attribute_names}
        node_data['sequence'] = seq
        node_data['direction'] = 1

        self.add_node(str(binode) + '_+', **node_data)

        node_data['sequence'] = sequence_complement(seq)
        node_data['direction'] = -1

        self.add_node(str(binode) + '_-', **node_data)

    def add_biedge(self, node1: str, node2: str, *, weight: int = 0):
        """
        Adds a biedge to the bidirected graph, comprising an edge and its complement.
        """
        if self.has_edge(node1, node2):
            raise ValueError(f'Attempted to add duplicate biedge: {node1}, {node2}')

        edge_data = {key: 0 for key in self.biedge_attribute_names}
        edge_data['weight'] = weight
        edge_data['index'] = self.number_of_biedges
        edge_data['is_representative'] = True
        self.add_edge(node1, node2, **edge_data)

        edge_data['is_representative'] = False
        self.add_edge(node_complement(node2), node_complement(node1), **edge_data)
        self.number_of_biedges += 1

    def representative_edge(self, edge: tuple):
        """
        Returns the edge associated with a biedge which is in the .gfa file.
        :param edge: one of the two edges in the biedge
        :return: the one which is in the .gfa file
        """
        return edge if self.edges[edge]['is_representative'] else edge_complement(edge)

    def reference_tree_edge(self, variant_edge):
        _, v = self.positive_variant_edge(variant_edge)
        if self.is_inversion(variant_edge):
            v = self.positive_node(v)
        w = self.parent_in_tree(v)
        return w, v

    def positive_variant_edge(self, edge: tuple):
        return edge if self.direction(edge[0]) == 1 or self.is_inversion(edge) else edge_complement(edge)

    def compute_edge_weights(self, walks: list[list]):
        """
        Computes the number of times that each edge is visited by some walk.
        """
        for walk in walks:
            for u, v in zip(walk[:-1], walk[1:]):
                self.edges[u, v]['weight'] += 1
                self.edges[node_complement(v), node_complement(u)]['weight'] += 1

    def annotate_branch_points(self) -> None:
        """
        Computes the branch point, i.e. the lowest common ancestor in the reference tree, of each variant biedge. The
        branch point is a positive orientation node. The edges (u,v), (u complement, v), (v,u), etc., all have the same
        branch point.
        """

        positive_direction_variants = [(self.positive_node(u), self.positive_node(v)) for u,v in self.variant_edges]

        # For each variant edge (u,v), yield the lowest common ancestor of u and v in the tree
        branch_point_tuples = nx.tree_all_pairs_lowest_common_ancestor(
            self.reference_tree,
            root=self.termini[0]+'_+',
            pairs=positive_direction_variants
        )
        branch_points = dict(branch_point_tuples)

        for edge in self.variant_edges:
            u, v = edge
            branch_point = branch_points[self.positive_node(u), self.positive_node(v)]
            if branch_point == v or branch_point == node_complement(u):
                assert self.reference_tree.in_degree(branch_point) == 1
                self.edges[edge]['is_back_edge'] = True
                self.edges[edge_complement(edge)]['is_back_edge'] = True
                branch_point = self.parent_in_tree(branch_point)
            self.edges[edge]['branch_point'] = branch_point
            self.edges[edge_complement(edge)]['branch_point'] = branch_point

    def walk_up_tree(self, ancestor, descendant, search_limit: int = inf) -> tuple[list, bool]:
        u = descendant
        result = []
        search_count = 0
        reach_limit = False
        while u != ancestor and search_count < search_limit:
            search_count += 1
            result.append(u)
            u = self.parent_in_tree(u)
        result.append(ancestor)
        if (search_count + 1) == search_limit:
            reach_limit = True
        return result, reach_limit

    def positive_node(self, node: str) -> str:
        return node if self.direction(node) == 1 else node_complement(node)

    def direction(self, node: str) -> int:
        return self.nodes[node]['direction']

    def get_variant_position(self, edge: tuple) -> int:
        u, v = edge
        if (self.is_snp(edge) or self.is_mnp(edge)) and self.nodes[v]['on_reference_path'] == 1:
            return int(self.nodes[u]['position']) + 1
        else:
            return int(self.nodes[u]['position'])

    def walk_sequence(self, walk: list[str]) -> str:
        seq = ''
        for node in walk:
            seq += self.nodes[node]['sequence']
        return seq

    def walk_with_variants(self, first: str, last: str, variant_edges: list, search_limit: int = inf) -> list:
        """
        Computes a walk from first node to last node, including all of the variant edges in the list.
        The walk proceeds from the 'end' of the first node to the 'start' of the last node, where the sequence associated
        with each node lies between its 'start' and 'end'. Thus, the first node is excluded from the walk, and the
        last node is usually excluded. However, the last node is included if the walk goes from + to -, because the
        'start' of the - node is the 'end' of its complementary + node.

        :param first: first node of the walk, either + or -
        :param last: last node
        :param variant_edges: variant edges in the order they are encountered, with the correct orientations
        :return: the nodes on the walk in order, excluding first and last
        """
        walk_direction = self.direction(first)
        source_nodes = [first] + [v for _, v in variant_edges]
        sink_nodes = [u for u, _ in variant_edges] + [last]
        result = []
        for source, sink in zip(source_nodes, sink_nodes):
            include_source = self.direction(source) == walk_direction
            source = source if self.direction(source) == walk_direction else node_complement(source)
            sink = sink if self.direction(source) == walk_direction else node_complement(source)
            pair = (source, self.positive_node(sink)) if walk_direction == 1 else \
                (self.positive_node(sink), self.positive_node(source))
            segment, _ = self.walk_up_tree(*pair, search_limit=search_limit)
            if not include_source:
                segment = segment[:-1] if walk_direction == 1 else segment[1:]
            segment = segment[::-1]
            result += segment if walk_direction == 1 else walk_complement(segment)

        # exclude first and sometimes last nodes
        if self.direction(last) == -1 and walk_direction == 1:
            return result[1:]
        return result[1:-1]

    def ref_alt_alleles(self, variant_edge: tuple,
                        walkup_limit: int = inf,
                        return_search_bool: bool=False,
                        ):
        """
        Computes the reference allele and alternative allele of the branch point for each variant edge.
        :param variant_edges: list of tuples (u,v).
        :return: dict mapping variant edges to tuples (ref, alt)
        """

        if self.is_in_tree(variant_edge):
            raise ValueError("Ref and alt alleles are only defined for variant edges")

        u, v = variant_edge
        branch_point = self.edges[u, v]['branch_point']
        first, last = (branch_point, v) if self.direction(u) == 1 else (u, branch_point)

        ref_path = self.walk_with_variants(first, last, [], search_limit=walkup_limit)
        alt_path = self.walk_with_variants(first, last, [variant_edge], search_limit=walkup_limit)

        ref_search_limit, alt_search_limit = False, False
        if alt_search_limit:
            alt_allele = '.'
        else:
            alt_allele = self.walk_sequence(alt_path)
        if ref_search_limit:
            ref_allele = '.'
        else:
            ref_allele = self.walk_sequence(ref_path)

        branch_sequence = self.nodes[branch_point]['sequence']
        if not branch_sequence:
            branch_sequence = 'N'
        last_letter_of_branch_point = branch_sequence[-1]

        if return_search_bool:
            return ref_allele, alt_allele, last_letter_of_branch_point, branch_point, ref_search_limit, alt_search_limit
        else:
            return ref_allele, alt_allele, last_letter_of_branch_point, branch_point

    def genotype_and_linear_coverage_by_sample(self, walks) -> tuple[dict, dict, list]:
        """
        Integrates the genotype of each sample from the genotype of each walk.
        :param walks: list of walks for each sample
        :return:
        """
        cr_dict_haplotype = dict()
        ca_dict_haplotype = dict()
        linear_coverages = []

        for walk in walks:
            cr_ca_dicts, linear_coverage = self.genotype(walk, return_linear_coverage=True)
            cr_dict_walk, ca_dict_walk = cr_ca_dicts[0], cr_ca_dicts[1]
            for edge, count in cr_dict_walk.items():
                if edge in cr_dict_haplotype:
                    cr_dict_haplotype[edge] += count
                else:
                    cr_dict_haplotype[edge] = count
            for edge, count in ca_dict_walk.items():
                if edge in ca_dict_haplotype:
                    ca_dict_haplotype[edge] += count
                else:
                    ca_dict_haplotype[edge] = count
            linear_coverages.append(linear_coverage)

        return cr_dict_haplotype, ca_dict_haplotype, linear_coverages


    def genotype(self, walk: list[str], return_linear_coverage: bool = False) -> Union[list[dict], tuple[dict, tuple]]:
        """
        Computes the number of time that a walk visits each variant edge.
        :param walk: list of nodes
        :param return_linear_coverage: if True, returns a tuple of the genotype dictionary and the min position/max right position of nodes on the walk
        :return: dictionary of edge-count pairs, optionally the linear coverage
        """

        # Append start and end nodes to walk
        start = [self.termini[0] + '_+' if self.direction(walk[0]) == 1 else self.termini[1] + '_-']
        end = [self.termini[1] + '_+' if self.direction(walk[-1]) == 1 else self.termini[0] + '_-']
        walk = start + walk + end

        if not hasattr(self, 'ref_edge_set'):
            self.ref_edge_set = {self.representative_edge(self.reference_tree_edge(var_edge)) for var_edge in self.sorted_variant_edges(exclude_terminus=True)}

        ref_edge_set = self.ref_edge_set

        count_ref = {}
        count_alt = {}
        min_pos = inf
        max_pos = -inf
        for e in zip(walk[:-1], walk[1:]):
            if not self.has_edge(*e):
                raise ValueError(f"Specified list contains edge {e} which is not present in the graph")

            if not self.edges[e]['is_representative']:
                e = edge_complement(e)

            if not self.is_terminal(e[0]) and min(*self.position(e)) < min_pos:
                min_pos = min(*self.position(e))
            if not self.is_terminal(e[1]) and max(*self.right_position(e)) > max_pos:
                max_pos = max(*self.right_position(e))

            if self.edges[e]['is_in_tree']:
                if e in ref_edge_set:
                    count_ref[e] = count_ref.get(e, 0) + 1
            else:
                count_alt[e] = count_alt.get(e, 0) + 1

        cr_ca_dicts = [count_ref, count_alt]

        return (cr_ca_dicts, (min_pos, max_pos)) if return_linear_coverage else cr_ca_dicts

    def count_edge_visits(self, genotype: dict) -> dict:
        """
        Computes the number of time that a walk visits every edge, given the number of times it visits every
        variant edge.
        :param genotype: for every variant edge that is visited, the number of visits
        :return: for every edge that is visited, the number of visits
        """

        sinks: list[str] = []
        sources: dict[str, int] = {}

        # Add variant edge endpoints as sources or sinks depending on their respective directions
        for variant_edge, visit_count in genotype.items():
            if variant_edge not in self.variant_edges:
                raise ValueError("geno dictionary contains a key which is not a variant edge")

            u, v = variant_edge
            new_sinks = []
            new_sources = []
            if self.direction(u) == 1:
                new_sinks.append(u)
            else:
                new_sources.append(node_complement(u))
            if self.direction(v) == 1:
                new_sources.append(v)
            else:
                new_sinks.append(node_complement(v))

            for w in new_sources:
                sources[w] = sources.get(w, 0) + 1

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
            while sources.get(current_node, 0) == 0:
                # If current_node is the root, it means that the input genotype was invalid
                if self.reference_tree.in_degree(current_node) == 0:
                    raise ValueError("The input genotype does not correspond to any valid walk")

                previous_node = current_node
                current_node = self.parent_in_tree(current_node)
                edge_representative = self.representative_edge((current_node, previous_node))
                edge_visits[edge_representative] = edge_visits.get(edge_representative, 0) + 1

            # when reaching the source, it is "used up"
            sources[current_node] -= 1

        return edge_visits

    def allele_length(self) -> dict:
        """
        Computes the allele length for each variant edge.
        :return: dictionary mapping variant edges to its own allele length, len(ref_allele) + len(alt_allele)
        """
        def _allele_length(variant_edge):
            ref_allele, alt_allele, _, _ = self.ref_alt_alleles(variant_edge)
            return len(ref_allele) + len(alt_allele)

        return {e: _allele_length(e) for e in self.sorted_variant_edges(exclude_terminus=True)}

    def allele_count(self) -> dict:
        """
        Computes alt and total allele counts, defined as the number of times that a walk visits a variant
        edge (u,v) and the corresponding branch point w.
        :return: dictionary mapping variant edges to (ref_count, alt_count) pairs
        """
        # TODO handles inversions correctly?
        return {e: (self.edges[self.reference_tree_edge(e)]['weight'], self.edges[e]['weight'])
                for e in self.sorted_variant_edges(exclude_terminus=True)}

    def compute_binode_positions(self):
        """
        Computes the position of each binode along the linear reference path, as well as the distance from the linear
        reference, in basepairs.
        """
        for node in self.reference_tree.nodes():
            self.nodes[node]['distance_from_reference'] = inf  # contigs not reachable from reference are at distance infinity

        current_position = 0
        for u in self.reference_path:
            current_position += len(self.nodes[u]['sequence'])
            self.nodes[u]['position'] = current_position
            self.nodes[u]['distance_from_reference'] = 0

        order = list(nx.topological_sort(self.reference_tree))
        for u in order[1:]: # skip the root
            predecessor = self.parent_in_tree(u)
            self.nodes[u]['position'] = np.maximum(self.nodes[u]['position'],
                                               self.nodes[predecessor]['position'])
            self.nodes[node_complement(u)]['position'] = self.nodes[u]['position']

            if self.nodes[u]['on_reference_path']:
                continue
            self.nodes[u]['distance_from_reference'] = (self.nodes[predecessor]['distance_from_reference'] +
                                                        len(self.nodes[u]['sequence']))
            self.nodes[node_complement(u)]['distance_from_reference'] = self.nodes[u]['distance_from_reference']

    def compute_binode_right_positions(self):
        positive_subgraph = self.subgraph([n for n, direction in self.nodes(data="direction") if direction == 1])
        positive_subgraph = self.edge_subgraph([edge for edge in positive_subgraph.edges() if not self.is_back_edge(edge)])
        order = nx.topological_sort(positive_subgraph)
        for u in reversed(list(order)):
            if self.on_reference_path(u):
                self.nodes[u]['right_position'] = self.nodes[u]['position']
                self.nodes[node_complement(u)]['right_position'] = self.nodes[node_complement(u)]['position']
                continue
            successor_positions = [self.right_position(v) for v in self.successors(u)]
            self.nodes[u]['right_position'] = np.min(successor_positions)
            self.nodes[node_complement(u)]['right_position'] = self.nodes[u]['right_position']

    def get_variants_at_interval(self, half_open_interval: tuple[int, int], exclude_root_edges=True) -> list:
        """
        Computes variant edges whose position intersects some interval, by iterating over all variant edges.
        :param half_open_interval: [start, end) of the interval
        :param exclude_root_edges: if True, edges beginning or ending at a terminus are excluded
        :return: list of edges
        """
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
            raise ValueError("Expected exactly 2 variants in the variant list")

        if any([self.is_back_edge(edge) for edge in variants]):
            return "Not triallelic"

        if self.is_inversion(variants[0]) != self.is_inversion(variants[1]):
            "Found an inversion and a non-inversion in the variant list"

        endpoint_nodes = [node + '_+' for node in endpoint_binodes]
        endpoint_nodes += [node + '_-' for node in endpoint_binodes]
        if not all([self.has_node(node) for node in endpoint_nodes]):
            raise ValueError("One or more of the endpoint nodes is not in the graph")

        # Detect which node of each binode demarcates the superbubble
        variant_branch_points = [self.edges[e]['branch_point'] for e in variants]
        assert not any([brach_point == 0 for brach_point in variant_branch_points]), "Branch point of variant not found"
        start_nodes = [node for node in endpoint_nodes if node in variant_branch_points]
        if len(start_nodes) == 0:
            return "Not triallelic"
        variant_end_points = []
        for u, v in variants:
            if self.direction(u) == 1 and self.direction(v) == 1:
                variant_end_points.append(v)
            elif self.direction(u) == -1 and self.direction(v) == -1:
                variant_end_points.append(node_complement(u))
            elif self.direction(u) == 1 and self.direction(v) == -1:
                variant_end_points.append(u)
                variant_end_points.append(node_complement(v))
            elif self.direction(u) == -1 and self.direction(v) == 1:
                variant_end_points.append(node_complement(u))
                variant_end_points.append(v)
            else:
                raise ValueError("Invalid direction for variant edge nodes.")
        # variant_end_points = [v if self.direction(u) == 1 else u for u, v in variants]
        end_nodes = [node for node in endpoint_nodes if node in variant_end_points]
        if len(end_nodes) == 0:
            return "Not triallelic"

        assert len(start_nodes) == 1 and len(end_nodes) == 1, \
            f"Found {len(start_nodes)} possible start nodes, and {len(end_nodes)} possible end nodes"

        start_node = start_nodes[0]
        end_node = end_nodes[0]

        # Get in-neighbors and out-neighbors excluding specific nodes
        start_degree = len({successor for successor in self.successors(start_node) if not self.is_terminal(successor)})
        end_degree = len({predecessor for predecessor in self.predecessors(end_node) if not self.is_terminal(predecessor)})

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


    def get_missing_variants(self, linear_coverages: list[tuple]) -> list:
        """
        Computes variant edges that are missing from a list of genotypes where each genotype represents one of the walks of a haplotype.
        :param genotypes: list of genotypes, where each genotype is a dictionary mapping variant edges to visit count"""

        # order walks and variants by position
        source_positions = np.sort([x[1] for x in linear_coverages] + [self.position('+_terminus_+')])
        sink_positions = np.sort([x[0] for x in linear_coverages] + [self.right_position('-_terminus_+')])
        sorted_variant_edges = self.sorted_variant_edges(exclude_terminus=True)
        sorted_variant_positions = [min(*self.position(e)) for e in sorted_variant_edges]

        result = []
        for source, sink in zip(source_positions, sink_positions):
            # indices of first and last variant edges u,v s.t. position of u in between source and sink
            first = np.searchsorted(sorted_variant_positions, source, side='left')
            last = np.searchsorted(sorted_variant_positions, sink, side='right')
            for i in range(first, last):
                if max(*self.right_position(sorted_variant_edges[i])) <= sink:
                    result.append(sorted_variant_edges[i])

        return result


    def match_sequence_up_tree(self, sequence: str, node: str) -> bool:
        node_sequence = self.nodes[node]['sequence']
        putative_match_length = min(len(sequence), len(node_sequence))
        if self.nodes[node]['sequence'][-putative_match_length:] != sequence[-putative_match_length:]:
            return False
        if putative_match_length == len(sequence):
            return True
        remaining_sequence = sequence[:putative_match_length]
        return self.match_sequence_up_tree(remaining_sequence, self.parent_in_tree(node))
        
    def match_sequence_down_tree(self, sequence: str, node: str) -> bool:
        node_sequence = self.nodes[node]['sequence']
        putative_match_length = min(len(sequence), len(node_sequence))
        if self.nodes[node]['sequence'][:putative_match_length] != sequence[:putative_match_length]:
            return False
        if putative_match_length == len(sequence):
            return True
        remaining_sequence = sequence[putative_match_length:]
        tree_successors = self.reference_tree.successors(node)
        return any([self.match_sequence_down_tree(remaining_sequence, successor) for successor in tree_successors])
        
    def annotate_repeat_motif(self, variant_edge: tuple[str, str]) -> Optional[str]:
        """
        Returns the repeat motif of a variant edge if it is a repeat. 
        Otherwise, returns None.
        """
        if self.is_inversion(variant_edge):
            return None
        variant_edge = self.positive_variant_edge(variant_edge)
        _, v = variant_edge
        ref_allele, alt_allele, _, branch_point = self.ref_alt_alleles(variant_edge)

        def get_repeat_motif(allele: str) -> Optional[str]:
            non_basepair_character = 'N'
            if any(letter == non_basepair_character for letter in allele):
                return None

            allele_length = len(allele)
            for repeat_length in range(1,allele_length+1):
                if allele_length % repeat_length != 0:
                    continue
                motif = allele[:repeat_length]
                if allele == motif * (allele_length // repeat_length):
                    break

            if self.match_sequence_up_tree(motif, branch_point):
                return motif
            if self.match_sequence_down_tree(motif, v):
                return motif
            return None

        
        if len(alt_allele) == 0:
            return get_repeat_motif(ref_allele)
        if len(ref_allele) == 0:
            return get_repeat_motif(alt_allele)
        














