import sys
from functools import lru_cache, cached_property
from linecache import cache
from math import inf
import networkx as nx
import numpy as np
from .utils import (
    read_gfa,
    read_gfa_line_by_line,
    node_complement,
    edge_complement,
    sequence_complement,
    walk_complement,
    group_walks_by_name,
    nearly_identical_alleles,
    _node_recover,
    merge_dicts,
    log_action
)
from .search_tree import assign_node_directions, max_weight_dfs_tree
import os
import time
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

    @cached_property
    def sorted_variant_edges(self) -> list[str]:
        sorted_vars = [edge for edge in self.variant_edges if not self.is_terminal(edge)]
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
            if not self.has_edge(*node_or_edge):
                raise ValueError("Graph does not have edge {node_or_edge}")
            return (self.nodes[node_or_edge[0]]['on_reference_path'] or self.nodes[node_or_edge[1]]['on_reference_path'])
        elif type(node_or_edge) is str:
            if not self.has_node(node_or_edge):
                raise ValueError("Graph does not have node {node_or_edge}")
            return self.nodes[node_or_edge]['on_reference_path']

    def parent_in_tree(self, node: str) -> str:
        if self.reference_tree.in_degree(node) == 0:
            return None
        return next(self.reference_tree.predecessors(node))

    def position(self, node_or_edge: Any) -> Any:
        if isinstance(node_or_edge, str):
            return self.nodes[node_or_edge]['position']
        elif isinstance(node_or_edge, tuple):
            return self.position(node_or_edge[0]), self.position(node_or_edge[1])
        else:
            raise TypeError

    def right_position(self, node_or_edge: Any) -> Any:
        if isinstance(node_or_edge, str):
            return self.nodes[node_or_edge]['right_position']
        elif isinstance(node_or_edge, tuple):
            return self.right_position(node_or_edge[0]), self.right_position(node_or_edge[1])
        else:
            raise TypeError

    @classmethod
    def from_gfa_line_by_line(cls,
                 gfa_file: str,
                 edgeinfo_file: str = None,
                 nodeinfo_file: str = None,
                 compressed: bool = False,
                 log_path: str = None
                 ):
        """
        Reads a .gfa file into a PangenomeGraph object.
        :param gfa_file: path to a file name ending in .gfa
        :param edgeinfo_file: path to a previously-computed .edgeinfo file, to avoid re-doing work
        :param nodeinfo_file: path to a previously-computed .nodeinfo file, to avoid re-doing work
        :param compressed: set to True in order to read a gzipped .gfa file
        """

        if not os.path.exists(gfa_file):
            raise FileNotFoundError(gfa_file)
        if edgeinfo_file:
            if not os.path.exists(edgeinfo_file):
                raise FileNotFoundError(edgeinfo_file)

        G = cls()

        gfa_basename = os.path.basename(gfa_file)
        walk_start_nodes = []
        walk_end_nodes = []

        start_time = time.time()
        print("Reading gfa file")
        for parts in read_gfa_line_by_line(gfa_file, compressed=compressed):
            if parts[0] == 'S':
                binode, sequence = parts[1], parts[2]
                G.add_binode(binode, sequence)
            elif parts[0] == 'L':
                biedge = parts[1]
                node1 = biedge[0] + '_' + biedge[2]
                node2 = biedge[1] + '_' + biedge[3]
                G.add_biedge(node1, node2)
            elif parts[0] == 'W':
                hit_reference = parts[1]
                # sample_name = parts[2]
                walk = parts[3]
                if hit_reference:
                    G.add_reference_path(walk)

                walk_start_nodes.append(walk[0])
                walk_end_nodes.append(walk[-1])
                if not edgeinfo_file:
                    G.compute_edge_weights([walk])

        if log_path:
            log_action(log_path, start_time, f"Reading gfa file: {gfa_basename}")

        print("Num of binodes:", (len(G.nodes) / 2))
        print("Num of biedges:", (len(G.edges) / 2))

        if nodeinfo_file:
            start_time = time.time()
            print("Reading nodeinfo file")
            G.read_nodeinfo(nodeinfo_file)
            if log_path:
                log_action(log_path, start_time, f"Reading nodeinfo file: {gfa_basename}")
        else:
            start_time = time.time()
            print("Assigning node directions")
            assign_node_directions(G, G.reference_path)
            if log_path:
                log_action(log_path, start_time, f"Assigning node directions: {gfa_basename}")

        G.add_terminal_nodes(walk_start_nodes=walk_start_nodes, walk_end_nodes=walk_end_nodes)
        if edgeinfo_file:
            start_time = time.time()
            print("Reading edgeinfo file")
            G.read_edgeinfo(edgeinfo_file)
            if log_path:
                log_action(log_path, start_time, f"Reading edgeinfo file: {gfa_basename}")
        else:
            start_time = time.time()
            print("Computing reference tree")
            G.compute_reference_tree()
            if log_path:
                log_action(log_path, start_time, f"Computing reference tree: {gfa_basename}")

            start_time = time.time()
            print("Computing branch points")
            G.annotate_branch_points()
            if log_path:
                log_action(log_path, start_time, f"Computing branch points: {gfa_basename}")

        if not nodeinfo_file:
            start_time = time.time()
            print("Computing positions")
            G.compute_binode_positions()
            G.compute_binode_right_positions()
            if log_path:
                log_action(log_path, start_time, f"Computing positions: {gfa_basename}")

        return G

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

    def identify_variant_type(self,
                              edge: tuple[str, str],
                              ref: str = None,
                              alt: str = None,
                              ) -> str:
        var_type = None
        if self.is_inversion(edge):
            if var_type is not None:
                raise KeyError(f'Variant, {edge}, has dual type.')
            var_type = 'inversion'
        if self.is_replacement(edge, ref=ref, alt=alt):
            if var_type is not None:
                raise KeyError(f'Variant, {edge}, has dual type.')
            var_type = 'replacement'
        if self.is_insertion(edge, ref=ref, alt=alt):
            if var_type is not None:
                raise KeyError(f'Variant, {edge}, has dual type.')
            var_type = 'insertion'
        if self.is_snp(edge, ref=ref, alt=alt):
            if var_type is not None:
                raise KeyError(f'Variant, {edge}, has dual type.')
            var_type = 'snp'
        if self.is_mnp(edge, ref=ref, alt=alt):
            if var_type is not None:
                raise KeyError(f'Variant, {edge}, has dual type.')
            var_type = 'mnp'
        if self.is_back_edge(edge):
            if var_type is not None:
                raise KeyError(f'Variant, {edge}, has dual type.')
            var_type = 'duplication'
        if self.is_forward_edge(edge):
            if var_type is not None:
                raise KeyError(f'Variant, {edge}, has dual type.')
            var_type = 'deletion'
        return var_type

    def is_inversion(self, edge: tuple[str, str]) -> bool:
        u, v = edge
        return self.direction(u) != self.direction(v)

    def is_in_tree(self, edge: tuple[str, str]) -> bool:
        return (edge in self.reference_tree or edge_complement(edge) in self.reference_tree)

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

    def is_insertion(self,
                     edge: tuple[str, str],
                     ref: str = None,
                     alt: str = None,
                     ) -> bool:
        if not self.is_crossing_edge(edge):
            return False
        if ref is None or alt is None:
            ref, _, _, _ = self.ref_alt_alleles(edge)
        return len(ref) == 0

    def is_replacement(self,
                       edge: tuple[str, str],
                       ref: str = None,
                       alt: str = None,
                       ) -> bool:
        if not self.is_crossing_edge(edge):
            return False
        if ref is None or alt is None:
            ref, alt, _, _ = self.ref_alt_alleles(edge)
        return len(ref) != 0 and len(alt) != 0 and len(ref) != len(alt)

    def is_snp(self,
               edge: tuple[str, str],
               ref: str = None,
               alt: str = None,
               ) -> bool:
        if not self.is_crossing_edge(edge):
            return False
        if ref is None or alt is None:
            ref, alt, _, _ = self.ref_alt_alleles(edge)
        return len(ref) == len(alt) == 1

    def is_mnp(self,
               edge: tuple[str, str],
               ref: str = None,
               alt: str = None,
               ) -> bool:
        if not self.is_crossing_edge(edge):
            return False
        if ref is None or alt is None:
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
                  gfa_path: str,
                  vcf_filename: str,
                  chr_name: str,
                  size_threshold: int = None,
                  check_degenerate: bool = False,
                  log_path: str = None
                  ) -> None:
        """
        Writes the variant call format (vcf) file.
        :param gfa_path: the .gfa file, from which walks are read
        :param vcf_filename: the output vcf file path
        :param chr_name: the chromosome name in the first column of output vcf file
        :param size_threshold: the truncation length of ref and alt sequence
        :param check_degenerate: whether to exclude variants whose ref and alt alleles are identical
        :return:
        """
        # 'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'sample1', 'sample2', ...
        meta_info = f'##fileformat=VCFv4.2\n'
        meta_info += f'##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        meta_info += f'##FORMAT=<ID=CR,Number=1,Type=String,Description="The REF count for haplotype.">\n'
        meta_info += f'##FORMAT=<ID=CA,Number=1,Type=String,Description="The ALT count for haplotype.">\n'
        meta_info += f'##INFO=<ID=OR,Number=1,Type=String,Description="Off-linear reference.">\n'
        meta_info += f'##INFO=<ID=VT,Number=1,Type=Integer,Description="Variant type.">\n'
        meta_info += f'##INFO=<ID=DR,Number=1,Type=Integer,Description="Distance of node u, v from reference.">\n'
        meta_info += f'##INFO=<ID=RC,Number=1,Type=String,Description="The REF allele count.">\n'
        meta_info += f'##INFO=<ID=AC,Number=1,Type=String,Description="The ALT allele count.">\n'
        meta_info += f'##INFO=<ID=AN,Number=1,Type=String,Description="The non-missing value count.">\n'
        meta_info += f'##INFO=<ID=PU,Number=1,Type=Integer,Description="Position of U (left node of variant edge)">\n'
        meta_info += f'##INFO=<ID=PV,Number=1,Type=Integer,Description="Position of V (right node of variant edge)">\n'
        meta_info += f'##INFO=<ID=TR_MOTIF,Number=1,Type=String,Description="Repeat motif">\n'
        meta_info += f'##INFO=<ID=NIA,Number=0,Type=Flag,Description="Nearly identical alleles">\n'
        meta_info += f'##contig=<ID={chr_name[3:]}>\n'

        gfa_basename = os.path.basename(gfa_path)
        allele_count_dict = self.allele_count()

        sample_cr_dict = defaultdict(dict)
        sample_ca_dict = defaultdict(dict)
        sample_lc_dict = defaultdict(list)

        sample_vcf_info_dict = dict()

        start_time = time.time()
        print("Computing genotype for haplotypes")
        pre_sample_name = None
        # Memory efficient way to read gfa data
        for parts in read_gfa_line_by_line(gfa_path):
            if parts[0] != 'W':
                continue
            haplotype_name = parts[2]
            sample_name = parts[2].split('_')[0]
            if pre_sample_name is not None and sample_name != pre_sample_name:
                # print("Sample:", pre_sample_name, sample_name)
                # print_current_memory_usage()
                assert sample_name not in sample_vcf_info_dict
                sample_missing_dict = {hap: set(self.get_missing_variants(lc)) for hap, lc in
                                       sample_lc_dict.items()}
                sample_info_list = self.get_sample_vcf_info(pre_sample_name,
                                                            sample_cr_dict,
                                                            sample_ca_dict,
                                                            sample_missing_dict)
                sample_vcf_info_dict[pre_sample_name] = sample_info_list

                sample_cr_dict.clear()
                sample_ca_dict.clear()
                sample_lc_dict.clear()

            walk = parts[3]
            cr_dict, ca_dict, linear_coverage = self.genotype(walk, return_linear_coverage=True)

            sample_cr_dict[haplotype_name] = merge_dicts([sample_cr_dict[haplotype_name], cr_dict])
            sample_ca_dict[haplotype_name] = merge_dicts([sample_ca_dict[haplotype_name], ca_dict])

            sample_lc_dict[haplotype_name].append(linear_coverage)

            pre_sample_name = sample_name

        assert sample_name not in sample_vcf_info_dict
        sample_missing_dict = {hap: set(self.get_missing_variants(lc)) for hap, lc in
                              sample_lc_dict.items()}
        sample_info_list = self.get_sample_vcf_info(sample_name,
                                                    sample_cr_dict,
                                                    sample_ca_dict,
                                                    sample_missing_dict
                                                    )
        sample_vcf_info_dict[sample_name] = sample_info_list

        sample_cr_dict.clear()
        sample_ca_dict.clear()
        sample_lc_dict.clear()

        if log_path:
            log_action(log_path, start_time, f"Reading gfa, computing genotype and missing variants for haplotypes: {gfa_basename}")

        # subject_ids = sorted({sample_name.split('_')[0] for sample_name in sample_walks_dict.keys() if not sample_name.startswith("GRCh38")})
        subject_ids = sorted(sample_vcf_info_dict.keys())

        header_names = list(self.vcf_attribute_names) + subject_ids

        start_time = time.time()
        print("Writing vcf file")
        with open(vcf_filename, 'w') as file:
            file.write(meta_info)
            file.write('#'+'\t'.join(header_names) + '\n')

            for idx, (u, v) in enumerate(tqdm(self.sorted_variant_edges)):
                representative_variant_edge = (u, v)
                # representative_ref_edge = self.representative_edge(self.reference_tree_edge(representative_variant_edge))

                if self.direction(u) == -1 and self.direction(v) == -1:
                    u, v = edge_complement((u, v))
                edge = (u, v)

                ref_allele, alt_allele, last_letter_of_branch_point, branch_point = self.ref_alt_alleles(edge)
                VT = self.identify_variant_type(edge, ref_allele, alt_allele)

                motif = self.annotate_repeat_motif(representative_variant_edge,
                                                   ref_allele=ref_allele,
                                                   alt_allele=alt_allele,
                                                   branch_point=branch_point)
                motif = '.' if motif is None else motif



                if check_degenerate:
                    if ref_allele == alt_allele:
                        continue

                prepend_letter_to_alleles = (len(ref_allele) == 0 or len(alt_allele) == 0)

                new_ref = '.'
                ref_allele_on_forward_reference_path = self.direction(edge[0]) == 1 and self.on_reference_path(edge)
                if not ref_allele_on_forward_reference_path:
                    new_ref = ref_allele
                    ref_allele = '.'
                    prepend_letter_to_alleles = False
                
                if prepend_letter_to_alleles:
                    ref_allele = last_letter_of_branch_point + ref_allele
                    alt_allele = last_letter_of_branch_point + alt_allele
                
                if size_threshold:
                    ref_allele = ref_allele[:size_threshold]
                    alt_allele = alt_allele[:size_threshold]

                edge_vcf_position = self.get_vcf_position(edge, prepend_letter_to_alleles)

                allele_data_list = []
                # 'CHROM' 0
                allele_data_list.append(chr_name)
                # 'POS' 1
                allele_data_list.append(str(edge_vcf_position))
                # 'ID' 2
                allele_data_list.append(''.join(tuple(map(lambda x: _node_recover(x), representative_variant_edge))))
                # 'REF' 3
                allele_data_list.append(ref_allele)
                # 'ALT' 4
                allele_data_list.append(alt_allele)
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
                    counts = sample_vcf_info_dict[sample_name][idx]
                    allele_data_list.append(counts)
                    for count in counts.split('|'):
                        if count != '.:.:.':
                            AN += 1

                RC = allele_count_dict[representative_variant_edge][0] if not self.is_inversion(edge) else '.'
                AC = allele_count_dict[representative_variant_edge][1]

                INFO = (f'OR={new_ref};'
                        f'VT={VT};'
                        f'DR={int(self.nodes[u]["distance_from_reference"])},{int(self.nodes[v]["distance_from_reference"])};'
                        f'RC={RC};'
                        f'AC={AC};'
                        f'AN={AN};'
                        f'PU={int(self.nodes[u]["position"])};'
                        f'PV={int(self.nodes[v]["position"])};'
                        f'TR_MOTIF={motif}')

                if nearly_identical_alleles(ref_allele, alt_allele):
                    INFO += ';NIA=1'
                else:
                    INFO += ';NIA=0'

                allele_data_list[7] = INFO

                file.write('\t'.join(allele_data_list) + '\n')
        if log_path:
            log_action(log_path, start_time, f"Writing vcf: {gfa_basename}")

    def get_sample_vcf_info(self,
                            sample_name,
                            sample_cr_dict,
                            sample_ca_dict,
                            sample_missing_dict
                            ):
        sample_vcf_info = []
        for u, v in self.sorted_variant_edges:
            representative_variant_edge = (u, v)
            representative_ref_edge = self.representative_edge(self.reference_tree_edge(representative_variant_edge))

            if self.direction(u) == -1 and self.direction(v) == -1:
                u, v = edge_complement((u, v))
            edge = (u, v)

            if sample_name.startswith(("CHM", "GRCh")):
                haplotype_name = sample_name + '_0'
                cr_0 = sample_cr_dict[haplotype_name].get(representative_ref_edge, 0) if not self.is_inversion(edge) else '.'
                ca_0 = sample_ca_dict[haplotype_name].get(representative_variant_edge, 0)

                counts = (f"{int(bool(ca_0))}:{cr_0}:{ca_0}"
                          if representative_variant_edge not in sample_missing_dict[haplotype_name] else '.:.:.')
                # counts = (f"{int(bool(ca_0))}:{cr_0}:{ca_0}"
                #           if cr_0 != 0 or ca_0 != 0 else '.:.:.')
            else:
                haplotype1_name = sample_name + '_1'
                haplotype2_name = sample_name + '_2'

                cr_1 = sample_cr_dict[haplotype1_name].get(representative_ref_edge, 0) if not self.is_inversion(edge) else '.'
                ca_1 = sample_ca_dict[haplotype1_name].get(representative_variant_edge, 0)

                cr_2 = sample_cr_dict[haplotype2_name].get(representative_ref_edge, 0) if not self.is_inversion(edge) else '.'
                ca_2 = sample_ca_dict[haplotype2_name].get(representative_variant_edge, 0)

                count_1 = (f"{int(bool(ca_1))}:{cr_1}:{ca_1}"
                           if representative_variant_edge not in sample_missing_dict[haplotype1_name] else '.:.:.')
                count_2 = (f"{int(bool(ca_2))}:{cr_2}:{ca_2}"
                           if representative_variant_edge not in sample_missing_dict[haplotype2_name] else '.:.:.')
                # count_1 = (f"{int(bool(ca_1))}:{cr_1}:{ca_1}"
                #            if cr_1 != 0 or ca_1 != 0 else '.:.:.')
                # count_2 = (f"{int(bool(ca_2))}:{cr_2}:{ca_2}"
                #            if cr_2 != 0 or ca_2 != 0 else '.:.:.')
                counts = f"{count_1}|{count_2}"
            sample_vcf_info.append(counts)
        return sample_vcf_info

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

    def compute_edge_weights(self, walks: list[list[str]]):
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

    def walk_up_tree(self, ancestor, descendant) -> list:
        u = descendant
        result = []
        while u != ancestor:
            result.append(u)
            u = self.parent_in_tree(u)
        result.append(ancestor)
        return result

    def positive_node(self, node: str) -> str:
        return node if self.direction(node) == 1 else node_complement(node)

    def direction(self, node: str) -> int:
        return self.nodes[node]['direction']

    def get_vcf_position(self,
                             edge: tuple,
                             prepend_letter_to_alleles: bool,
                             ) -> int:
        """The VCF position of a variant is offset by 1 compared with the ordinary position, except
        variants that by convention have the last letter of the branch point prepended to their ref 
        and their alt allele."""
        u, _ = self.positive_variant_edge(edge)
        return self.position(u) + 1 - int(prepend_letter_to_alleles)

    def walk_sequence(self, walk: list[str]) -> str:
        seq = ''
        for node in walk:
            seq += self.nodes[node]['sequence']
        return seq

    def walk_with_variants(self, first: str, last: str, variant_edges: list) -> list:
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
            segment = self.walk_up_tree(*pair)
            if not include_source:
                segment = segment[:-1] if walk_direction == 1 else segment[1:]
            segment = segment[::-1]
            result += segment if walk_direction == 1 else walk_complement(segment)

        # exclude first and sometimes last nodes
        if self.direction(last) == -1 and walk_direction == 1:
            return result[1:]
        return result[1:-1]

    def ref_alt_alleles(self,
                        variant_edge: tuple
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

        ref_path = self.walk_with_variants(first, last, [])
        alt_path = self.walk_with_variants(first, last, [variant_edge])

        alt_allele = self.walk_sequence(alt_path)
        ref_allele = self.walk_sequence(ref_path)

        branch_sequence = self.nodes[branch_point]['sequence']
        if not branch_sequence:
            branch_sequence = 'N'
        last_letter_of_branch_point = branch_sequence[-1]

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


    def genotype(self, walk: list[str], return_linear_coverage: bool = False):
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
            self.ref_edge_set = {self.representative_edge(self.reference_tree_edge(var_edge)) for var_edge in self.sorted_variant_edges}

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

        return (count_ref, count_alt, (min_pos, max_pos)) if return_linear_coverage else (count_ref, count_alt)

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

        return {e: _allele_length(e) for e in self.sorted_variant_edges}

    def allele_count(self) -> dict:
        """
        Computes alt and total allele counts, defined as the number of times that a walk visits a variant
        edge (u,v) and the corresponding branch point w.
        :return: dictionary mapping variant edges to (ref_count, alt_count) pairs
        """
        # TODO handles inversions correctly?
        return {e: (self.edges[self.reference_tree_edge(e)]['weight'], self.edges[e]['weight'])
                for e in self.sorted_variant_edges}

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
        """Computes the right position of each binode, defined as the minimum position of its successors in the 
        positive subgraph minus back edges."""
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
            return "not_triallelic"

        if self.is_inversion(variants[0]) != self.is_inversion(variants[1]):
           raise ValueError("Found an inversion and a non-inversion in the variant list")

        endpoint_nodes = [node + '_+' for node in endpoint_binodes]
        endpoint_nodes += [node + '_-' for node in endpoint_binodes]
        if not all([self.has_node(node) for node in endpoint_nodes]):
            raise ValueError("One or more of the endpoint nodes is not in the graph")

        # Detect which node of each binode demarcates the superbubble
        variant_branch_points = [self.edges[e]['branch_point'] for e in variants]
        assert not any([brach_point == 0 for brach_point in variant_branch_points]), "Branch point of variant not found"
        start_nodes = [node for node in endpoint_nodes if node in variant_branch_points]
        if len(start_nodes) == 0:
            return "not_triallelic"
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

        end_nodes = [node for node in endpoint_nodes if node in variant_end_points]
        if len(end_nodes) == 0:
            return "not_triallelic"

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
        Computes variant edges that are missing from a haplotype.
        :param linear_coverages: minimum and maximum positions of each walk in a haplotype."""

        # order walks and variants by position
        source_positions = np.sort([x[1] for x in linear_coverages] + [self.position('+_terminus_+')])
        sink_positions = np.sort([x[0] for x in linear_coverages] + [self.right_position('-_terminus_+')])
        sorted_variant_edges = self.sorted_variant_edges
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


    def _match_sequence_up_tree(self, sequence: str, node: str) -> bool:
        """Returns True if `sequence` is a suffix of the unique sequence beginning at a terminus and
        ending at `node` within the tree.
        """
        node_sequence = self.nodes[node]['sequence']
        putative_match_length = min(len(sequence), len(node_sequence))
        if self.nodes[node]['sequence'][-putative_match_length:] != sequence[-putative_match_length:]:
            return False
        if putative_match_length == len(sequence):
            return True
        remaining_sequence = sequence[:putative_match_length]
        return self._match_sequence_up_tree(remaining_sequence, self.parent_in_tree(node))
        
    def _match_sequence_down_tree(self, sequence: str, node: str) -> bool:
        """Returns True if `sequence' is a prefix of some sequence beginning at `node` in
        thre tree."""
        node_sequence = self.nodes[node]['sequence']
        putative_match_length = min(len(sequence), len(node_sequence))
        if self.nodes[node]['sequence'][:putative_match_length] != sequence[:putative_match_length]:
            return False
        if putative_match_length == len(sequence):
            return True
        remaining_sequence = sequence[putative_match_length:]
        tree_successors = self.reference_tree.successors(node)
        return any([self._match_sequence_down_tree(remaining_sequence, successor) \
            for successor in tree_successors])
        
    def _match_sequence_down_graph(self, sequence: str, node: str, end_node: str) -> bool:
        """Returns True if `sequence' is a prefix of some sequence beginning at `node` in
        the graph."""
        node_sequence = self.nodes[node]['sequence']
        putative_match_length = min(len(sequence), len(node_sequence))
        if self.nodes[node]['sequence'][:putative_match_length] != sequence[:putative_match_length]:
            return False
        if putative_match_length == len(sequence):
            return (node == end_node and len(sequence) == len(node_sequence))
        remaining_sequence = sequence[putative_match_length:]
        tree_successors = self.successors(node)
        return any([self._match_sequence_down_graph(remaining_sequence, successor, end_node) \
            for successor in tree_successors])
        
    def annotate_repeat_motif(self,
                              variant_edge: tuple[str, str],
                              ref_allele: str = None,
                              alt_allele: str = None,
                              branch_point: str = None) -> Optional[str]:
        """
        Returns the repeat motif of a variant edge if it is a repeat. 
        Otherwise, returns None.
        """
        import sys
        sys.setrecursionlimit(len(self.reference_tree.nodes))
        if self.is_inversion(variant_edge):
            return None
        variant_edge = self.positive_variant_edge(variant_edge)
        _, v = variant_edge
        if ref_allele is None or alt_allele is None or branch_point is None:
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

            if self._match_sequence_up_tree(motif, branch_point):
                return motif
            if not self.is_back_edge(variant_edge):
                if self._match_sequence_down_tree(motif, v):
                    return motif
            else:
                return motif
            return None

        
        if len(alt_allele) == 0:
            return get_repeat_motif(ref_allele)
        if len(ref_allele) == 0:
            return get_repeat_motif(alt_allele)

    def missing_inversion_allele(self, variant_edge: tuple[str, str], minimum_alt_length: int=10) -> Optional[str]:
        """Checks if the alt allele of a variant edge matches the reverse complement of some other path from branch 
        point to v, whether that be the reference allele or a different alt path. Returns the matching allele
        or None."""
        if not self.is_crossing_edge(variant_edge):
            return None
        variant_edge = self.positive_variant_edge(variant_edge)
        ref, alt, _, branch_point = self.ref_alt_alleles(variant_edge)
        if len(alt) < minimum_alt_length:
            return None

        end_node = variant_edge[1]
        alt_complement = self.nodes[branch_point]['sequence'] \
            + sequence_complement(alt) + self.nodes[end_node]['sequence']

        if self._match_sequence_down_graph(alt_complement, branch_point, end_node):
            return alt

        return None
        
    def simplify_subgraph(self,
                          pos_range: tuple = None,
                          endpoints: set[set] = None,
                          minimum_allele_length: int = 1000) -> nx.DiGraph:
        """Returns a simplified minor of the underlying directed graph, with variant
        edges remvoed if their combined allele length is less than the minimum, with
        'tips' removed, and with 'paths' contracted."""

        subgraph = self.delete_small_variants(pos_range, minimum_allele_length)
        assert not nx.is_empty(subgraph), "Invalid position range: Empty subgraph."
        if endpoints is None:
            node_in_degrees = {node: degree for node, degree in subgraph.in_degree()}
            node_out_degrees = {node: degree for node, degree in subgraph.out_degree()}
            endpoints = {node for node in subgraph.nodes
                         if self.nodes[node]['on_reference_path'] == 1
                         and (node_in_degrees[node] == 0 or node_out_degrees[node] == 0)}
            assert len(endpoints) == 4
        self.delete_tips(subgraph, endpoints)
        self.contract_paths(subgraph, endpoints)

        return subgraph

    def delete_small_variants(self, pos_range: tuple = None, minimum_allele_length: int = 5) -> nx.DiGraph:
        """Returns an edge subgraph of the underlying directed graph with variant edges removed
        if their combined allele length is less than the minimum."""
        if pos_range == None:
            simplified_graph = nx.DiGraph(self)
            variant_edge_set = self.variant_edges
        else:
            nodes_to_include = [node for node, data in self.nodes(data=True)
                                if (data['position'] >= pos_range[0] and data['position'] < pos_range[1]) or
                                (data['position'] == pos_range[1] and data['on_reference_path'] == 1)]
            simplified_graph = nx.DiGraph(self.subgraph(nodes_to_include))
            variant_edge_set = [variant for variant in self.variant_edges
                                 if min(self.position(variant)) >= pos_range[0] and max(self.position(variant)) <= pos_range[1]]

        small_variant_edges = []
        for variant_edge in variant_edge_set:
            ref, alt, _, _ = self.ref_alt_alleles(variant_edge)
            if len(ref) + len(alt) < minimum_allele_length:
                small_variant_edges.append(variant_edge)
                small_variant_edges.append(edge_complement(variant_edge))

        simplified_graph.remove_edges_from(small_variant_edges)

        return simplified_graph

    @staticmethod
    def delete_tips(simplified_graph: nx.DiGraph, end_points: set[str]) -> None:
        node_in_degrees = {node: degree for node, degree in simplified_graph.in_degree()}
        node_out_degrees = {node: degree for node, degree in simplified_graph.out_degree()}
        tips = []
        for node, in_degree in node_in_degrees.items():
            out_degree = node_out_degrees[node]
            if in_degree == 1 and out_degree == 0 and node not in end_points:
                tips.append(node)

        nodes_to_delete = []
        for tip in tips:
            node = tip
            while node_in_degrees[node] == 1 and node_out_degrees[node] == 0:
                parent = next(simplified_graph.predecessors(node))
                node_out_degrees[parent] -= 1
                nodes_to_delete.append(node)
                nodes_to_delete.append(node_complement(node))
                node = parent

        simplified_graph.remove_nodes_from(nodes_to_delete)

    
    @staticmethod
    def contract_paths(simplified_graph: nx.DiGraph, end_points: set[str]) -> dict[str, str]:
        """
        Contracts paths in a graph by combining nodes.
        :param simplified_graph: A NetworkX directed graph with tips deleted and small variants removed
        :param end_points: A set of end points in the graph, which won't be combined with other nodes
        :return: A dictionary mapping each node that was combined with another node to the node
        with which it was combined
        """

        def combine_nodes(parent: str, child: str) -> None:
            neighbors = list(simplified_graph.successors(child))
            sequence = simplified_graph.nodes[child]['sequence']
            simplified_graph.remove_node(child)
            simplified_graph.remove_node(node_complement(child))
            for neighbor in neighbors:
                simplified_graph.add_edge(parent, neighbor)
                simplified_graph.add_edge(*edge_complement((parent, neighbor)))
            simplified_graph.nodes[parent]['sequence'] += sequence
            simplified_graph.nodes[node_complement(parent)]['sequence'] = sequence + \
                simplified_graph.nodes[node_complement(parent)]['sequence']
            

        node_in_degrees = {node: degree for node, degree in simplified_graph.in_degree()}
        node_out_degrees = {node: degree for node, degree in simplified_graph.out_degree()}

        starting_points = list(end_points)
        for node, out_degree in node_out_degrees.items():
            in_degree = node_in_degrees[node]
            if simplified_graph.nodes[node]['direction'] == 1 and (out_degree > 1 or in_degree > 1):
                starting_points.append(node)

        node_mapping: dict[str, str] = {}
        for start in starting_points:
            for node in list(simplified_graph.successors(start)):
                parent = start
                combined_nodes = []
                while node_in_degrees[node] == 1 and node_out_degrees[node] == 1:
                    if node_out_degrees[parent] == 1:
                        combine_nodes(parent, node) # TODO replace with combine_path
                        combined_nodes.append(node)
                        node = next(simplified_graph.successors(parent))
                    else:
                        parent = node
                        node = next(simplified_graph.successors(node))

                for combined_node in combined_nodes:
                    node_mapping[combined_node] = parent

        return node_mapping


        
            
        


        














