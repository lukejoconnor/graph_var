import os.path

from graph import PangenomeGraph
import re
import numpy as np
from tqdm import tqdm
from intervaltree import Interval, IntervalTree
from typing import List, Tuple, Union, Dict, Set, Optional
from collections import defaultdict
from utils import _node_convert, load_graph_from_pkl, save_graph_to_pkl
import pandas as pd
import ast

def get_node_id_symbol(node: str) -> Tuple[str, str]:
    node_split = node.split('_')
    if len(node_split) == 2:
        node_id, symbol = node_split
    elif len(node_split) == 3:  # for terminals
        direction, node_id, symbol = node_split
        node_id = direction + "_" + node_id
    else:
        raise ValueError(f"Invalid node ID: {node}")

    return node_id, symbol

def write_dfs_tree_to_gfa(G: PangenomeGraph, filename: str):
    with open(filename, 'wb') as gfa_file:
        for node in G.reference_tree.nodes(data=False):
            node_id, symbol = get_node_id_symbol(node)
            gfa_file.write(f'S\t{node_id}\t{G.nodes[node].get("sequence", "*")}\n'.encode())

        for edge in G.reference_tree.edges(data=False):
            u, v = edge
            u_id, u_symbol = get_node_id_symbol(u)
            v_id, v_symbol = get_node_id_symbol(v)
            gfa_file.write(f'L\t{u_id}\t{u_symbol}\t{v_id}\t{v_symbol}\t0M\n'.encode())

def write_nodes_to_txt(nodes: List[str], filename: str):
    with open(filename, 'wb') as txt_file:
        for node in nodes:
            txt_file.write(node.encode())
            txt_file.write('\n'.encode())

def extract_bubble_ids(node_string: str) -> Tuple:
    # Split the string and keep '>' and '<' symbols
    node_ids = re.findall(r'[><]\d+', node_string)

    # Convert the list of strings to a tuple
    node_tuple = tuple(map(lambda x: _node_convert(x), node_ids))

    assert len(node_tuple) == 2, f"Invalid bubble ID: {node_string}"

    return node_tuple

def extract_nodes_in_bubble(node_string: str) -> List[str]:
    # Split the string and keep '>' and '<' symbols
    node_ids = re.findall(r'[><]\d+', node_string)

    # Convert the list of strings to a tuple
    node_list = list(map(lambda x: _node_convert(x), node_ids))

    return node_list

def extract_node_bubble_partition_from_vcf(vcf_path: str) -> Tuple[Dict, Dict]:
    bubbles = dict()
    node_partition = defaultdict(set)
    with open(vcf_path, 'r') as vcf_file:
        for line in vcf_file:
            if line.startswith('##'):
                continue

            if line.startswith('#'):
                continue

            parts = line.strip().split('\t')
            ID = parts[2]
            INFO = parts[7]
            INFO_list = INFO.split(';')
            data_dict = {x.split('=')[0]: x.split('=')[1] for x in INFO_list}
            AT = data_dict['AT']
            bubble_id_tuple = extract_bubble_ids(ID)
            nodes_list = extract_nodes_in_bubble(AT)

            bubbles[bubble_id_tuple] = data_dict

            for node in nodes_list:
                node_partition[node].add(bubble_id_tuple)

    return bubbles, node_partition

def find_bubbles_from_vcf(G: PangenomeGraph, vcf_path: str) -> Dict:
    bubbles = {}
    with open(vcf_path, 'r') as vcf_file:
        for line in vcf_file:
            if line.startswith('##'):
                continue

            if line.startswith('#'):
                continue

            parts = line.strip().split('\t')
            ID = parts[2]
            bubble_id_tuple = extract_bubble_ids(ID)
            bubbles[bubble_id_tuple] = (G.nodes[bubble_id_tuple[0]]['position'], G.nodes[bubble_id_tuple[1]]['position'])

    sorted_dict = dict(sorted(bubbles.items(), key=lambda item: item[1][0]))

    return sorted_dict


def find_node_partition(G: PangenomeGraph, bubble_dict: Dict) -> Dict:
    # Small epsilon value to shift the start of intervals
    epsilon = 1e-1

    node_partition = {}
    position_sorted_nodes = sorted(G.nodes(data=True), key=lambda x: x[1]['position'])
    # Create the interval tree
    close_position_interval_tree = IntervalTree(Interval(min(v), max(v)+epsilon, k) for k, v in bubble_dict.items())

    for node, node_attr in position_sorted_nodes:
        node_partition[node] = set([interval.data for interval in close_position_interval_tree[node_attr['position']]])

    return node_partition


def get_variants_for_bubbles(G: PangenomeGraph,
                             node_partition: Dict[str, Set[Tuple[str, str]]],
                             inversion_only: bool = False,
                             ) -> Tuple[Dict, Dict]:
    bubble_within_variants = defaultdict(set)
    bubble_crossing_variants = defaultdict(set)

    for edge in G.variant_edges:
        u, v = edge

        if inversion_only:
            if not G.is_inversion(edge):
                continue

        bubbles_u = node_partition[u]
        bubbles_v = node_partition[v]

        bubble_set_intersection = bubbles_u.intersection(bubbles_v)
        bubble_set_complementary_u = bubbles_u - bubble_set_intersection
        bubble_set_complementary_v = bubbles_v - bubble_set_intersection

        if len(bubble_set_intersection) > 0:
            for bubble in bubble_set_intersection:
                bubble_within_variants[bubble].add(edge)
        else:
            for bubble in bubble_set_complementary_u:
                bubble_crossing_variants[bubble].add(edge)

            for bubble in bubble_set_complementary_v:
                bubble_crossing_variants[bubble].add(edge)

    return bubble_within_variants, bubble_crossing_variants

def variant_edges_summary(G: PangenomeGraph, var_list: list) -> Dict:
    summary_dict = dict()
    for edge in sorted(list(var_list)):
        if G.is_inversion(edge):
            summary_dict['inversions'] = summary_dict.get('inversions', 0) + 1
        if G.is_snp(edge):
            summary_dict['snps'] = summary_dict.get('snps', 0) + 1
        if G.is_mnp(edge):
            summary_dict['mnps'] = summary_dict.get('mnps', 0) + 1
        if G.is_crossing_edge(edge):
            summary_dict['crossing_edges'] = summary_dict.get('crossing_edges', 0) + 1
        if G.is_back_edge(edge):
            summary_dict['back_edges'] = summary_dict.get('back_edges', 0) + 1
        if G.is_forward_edge(edge):
            summary_dict['forward_edges'] = summary_dict.get('forward_edges', 0) + 1
    summary_dict['total'] = len(var_list)
    return summary_dict

def write_bubble_summary_result(gfa_path: str,
                                vcf_path: str,
                                save_pangenome_graph: Optional[bool] = False,
                                G_pkl_path: Optional[str] = None,
                                reference_path_index: Optional[int] = None,
                                output_dir: Optional[str] = './',
                                method: Optional[str] = "AT",
                                bubble_dict: Optional[Dict] = None,
                                node_partition: Optional[Dict] = None,
                                compressed: Optional[bool] = False):
    if not os.path.isfile(gfa_path):
        raise FileNotFoundError(f"GFA file not found: {gfa_path}")

    if not os.path.isfile(vcf_path):
        raise FileNotFoundError(f"VCF file not found: {vcf_path}")

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if G_pkl_path and os.path.isfile(G_pkl_path):
        print(f"Loading graph from {G_pkl_path}...")
        G, walks, walk_sample_names = load_graph_from_pkl(G_pkl_path, compressed=compressed)
        chr_name = os.path.basename(G_pkl_path).split('.')[0]
    else:
        print("Constructing graph from GFA...")
        G, walks, walk_sample_names = PangenomeGraph.from_gfa(gfa_path,
                                                              reference_path_index=reference_path_index,
                                                              return_walks=True,
                                                              compressed=False)
        chr_name = os.path.basename(gfa_path).split('.')[0]
        if save_pangenome_graph:
            if compressed:
                save_path = os.path.join(output_dir, chr_name+'.pkl.gz')
            else:
                save_path = os.path.join(output_dir, chr_name+'.pkl')
            save_graph_to_pkl(G, walks, walk_sample_names, save_path, compressed=compressed)

    print("Assigning node to bubbles...")

    if method == "AT":
        if node_partition is None or bubble_dict is None:
            bubble_dict, node_partition = extract_node_bubble_partition_from_vcf(vcf_path)
        method_suffix = "_AT"
    elif method == "Position":
        if node_partition is None or bubble_dict is None:
            bubble_dict = find_bubbles_from_vcf(G, vcf_path)
            node_partition = find_node_partition(G, bubble_dict)
        method_suffix = "_Position"
    else:
        raise ValueError(f"Invalid method: {method}")


    bubble_var_count_path = f"bubble_variant_counts_{chr_name}{method_suffix}.tsv"
    dfs_tree_path = f"dfs_tree_{chr_name}.gfa"

    print("Writing DFS tree to GFA...")
    write_dfs_tree_to_gfa(G, os.path.join(output_dir, dfs_tree_path))

    print("Writing bubble summary to CSV...")
    var_dict_within, var_dict_missing = get_variants_for_bubbles(G, node_partition)

    bubble_list = list(bubble_dict.keys())

    var_with = [var_dict_within.get(key, {}) for key in bubble_list]
    var_with_summary = [variant_edges_summary(G, var_dict_within.get(key, [])) for key in bubble_list]
    var_cross = [var_dict_missing.get(key, {}) for key in bubble_list]
    var_cross_summary = [variant_edges_summary(G, var_dict_missing.get(key, [])) for key in bubble_list]

    length_with = [len(var_dict_within.get(key, {})) for key in bubble_list]
    length_cross = [len(var_dict_missing.get(key, {})) for key in bubble_list]
    length_total = [length_with[i] + length_cross[i] for i in range(len(bubble_list))]

    if method == "AT":
        AC_sum = [sum(ast.literal_eval(f"[{bubble_dict[x]['AC']}]")) for x in bubble_list]
    else:
        AC_sum = ['.'] * len(bubble_list)

    count_summary_df = pd.DataFrame({
         "Bubble_ids": bubble_list,
         "Total_count": length_total,
         "AC_sum": AC_sum,
         "Within_count": length_with,
         "Crossing_count": length_cross,
         "Within": var_with,
         "Within_summary": var_with_summary,
         "Crossing": var_cross,
         "Crossing_summary": var_cross_summary,
         }
    )

    count_summary_df.to_csv(os.path.join(output_dir, bubble_var_count_path), sep='\t')




"""
def construct_intervaltree_variant_edges(G: PangenomeGraph) -> IntervalTree:
    interval_tree = IntervalTree()
    for edge in G.variant_edges:
        u, v = edge
        positions = [G.nodes[u]['position'], G.nodes[v]['position']]
        x, y = min(positions), max(positions)

        interval_tree.add(Interval(x, y, edge))

    return interval_tree

def get_variants_from_intervaltree(half_open_interval: Union[Tuple[int, int], List[Tuple[int, int]]],
                                   vertice_bubble_dict: dict,
                                   exclude_root_edges=True) -> dict:
    result = {}

    if isinstance(half_open_interval, tuple):
        half_open_interval = [half_open_interval]
    elif isinstance(half_open_interval, list):
        pass
    else:
        raise ValueError("Invalid input type for half_open_interval")

    for start, end in half_open_interval:
        positions = [self.nodes[u]['position'], self.nodes[v]['position']]
        x, y = min(positions), max(positions)
        if x < end and y >= start:
            result.append(edge)

    if exclude_root_edges:
        result = [(u, v) for u, v in result if u != '+_terminus_+' and v != '+_terminus_-']

    return result
"""
