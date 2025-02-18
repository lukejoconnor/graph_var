import os.path

from .graph import PangenomeGraph
from .utils import _node_convert, load_graph_from_pkl, save_graph_to_pkl
import re
import numpy as np
from tqdm import tqdm
from intervaltree import Interval, IntervalTree
from typing import List, Tuple, Union, Dict, Set, Optional
from collections import defaultdict
import pandas as pd
import ast

# Utility functions ----------------------------------------------------------------------------------------------------
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

def get_variants_for_bubbles(G: PangenomeGraph,
                             node_partition: Dict[str, Set[Tuple[str, str]]]
                             ) -> Tuple[Dict, Dict]:
    bubble_within_variants = defaultdict(set)
    bubble_crossing_variants = defaultdict(set)

    for edge in sorted(list(G.variant_edges)):
        u, v = edge

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

def get_variants_for_bubbles_from_vcf(
                             vcf_path: str,
                             node_partition: Dict[str, Set[Tuple[str, str]]]
                             ) -> Tuple[Dict, Dict]:
    bubble_within_variants = defaultdict(set)
    bubble_crossing_variants = defaultdict(set)

    vcf_df = read_vcf_to_dataframe(vcf_path)

    for i in range(len(vcf_df)):
        edge = eval(vcf_df.loc[i, 'ID'])
        u, v = edge

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

def get_interval_tree_from_bed(bed_file: str, chr_name: str) -> IntervalTree:
    interval_tree = IntervalTree()
    with open(bed_file, 'r') as bed:
        for line in bed:
            parts = line.strip().split('\t')
            if parts[0] != chr_name:
                continue
            start = int(parts[1])
            end = int(parts[2])
            interval_tree.add(Interval(start, end))
    return interval_tree

def read_vcf_to_dataframe(vcf_path: str, return_meta_info: bool = False):
    records = []
    meta_info = []

    # Loop through each record in the VCF file
    with open(vcf_path, 'r') as vcf_file:
        header = None
        for line in vcf_file:
            if line.startswith('##'):
                meta_info.append(line)
                continue

            if line.startswith('#'):
                header = line.strip().split('\t')
                continue

            parts = line.strip().split('\t')
            assert len(parts) == len(header), f"Invalid VCF format: {line}"
            rec_dict = {header[i]: parts[i] for i in range(len(header))}
            records.append(rec_dict)

    # Convert to DataFrame
    df = pd.DataFrame(records)

    if return_meta_info:
        return df, meta_info
    else:
        return df

def get_info_dict(info_str: str) -> Dict:
    info_list = info_str.split(';')
    info_dict = {x.split('=')[0]: x.split('=')[1] for x in info_list if x != ''}
    return info_dict

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

def write_node_sequence_to_csv(G: PangenomeGraph, filename: str):
    nodes = sorted([node[:-2] for node in list(G.nodes)])
    sequences = [G.nodes[G.positive_node(node+'_+')]['sequence'] for node in nodes]
    df = pd.DataFrame({'NodeID': nodes, 'Sequence': sequences})
    df.to_csv(filename, index=False)

def write_nodes_to_txt(nodes: List[str], filename: str):
    with open(filename, 'wb') as txt_file:
        for node in nodes:
            txt_file.write(node.encode())
            txt_file.write('\n'.encode())

def extract_bubble_ids(node_string: str, symbol=False) -> Tuple:
    # Split the string and keep '>' and '<' symbols
    node_ids = re.findall(r'[><]\d+', node_string)

    # Convert the list of strings to a tuple
    if symbol:
        node_tuple = tuple(map(lambda x: _node_convert(x), node_ids))
    else:
        node_tuple = tuple(map(lambda x: x[1:], node_ids))

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
            if line.startswith('#') or line.startswith('##'):
                continue

            parts = line.strip().split('\t')
            POS = parts[1]
            ID = parts[2]
            ref, alt = parts[3], parts[4]
            INFO = parts[7]
            INFO_list = INFO.split(';')
            data_dict = {x.split('=')[0]: x.split('=')[1] for x in INFO_list}
            data_dict['POS'] = POS
            data_dict['REF'] = ref
            data_dict['ALT'] = alt
            AT = data_dict['AT']
            bubble_id_tuple = extract_bubble_ids(ID)
            nodes_list = extract_nodes_in_bubble(AT)

            bubbles[bubble_id_tuple] = data_dict

            for node in nodes_list:
                node_partition[node].add(bubble_id_tuple)

    return bubbles, node_partition

def extract_node_bubble_partition_from_snarl(snarl_path: str) -> Tuple[Dict, Dict]:
    bubbles = dict()
    node_partition = defaultdict(set)
    with open(snarl_path, 'r') as table_file:
        for line in table_file:
            if line.startswith('#') or line.startswith('##'):
                continue

            parts = line.strip().split('\t')
            ID = parts[0]
            nodes = parts[3]
            bubble_id_tuple = extract_bubble_ids(ID)
            nodes_list = (list(map(lambda x: x + '_+', bubble_id_tuple)) +
                          list(map(lambda x: x + '_-', bubble_id_tuple)) +
                          list(map(lambda x: x + '_+', nodes.split(','))) +
                          list(map(lambda x: x + '_-', nodes.split(','))))

            bubbles[bubble_id_tuple] = {'Level': parts[1], 'Parent': parts[2], 'Content': parts[3]}

            for node in nodes_list:
                node_partition[node].add(bubble_id_tuple)

    return bubbles, node_partition

def find_variants_outside_interaltree(interval_tree: IntervalTree,
                                      G: PangenomeGraph,
                                      exclude_terminus: bool = False) -> list:
    egde_outside_region = []
    for edge in sorted(list(G.variant_edges)):
        u, v = edge

        if exclude_terminus:
            nodes_to_exclude = {'+_terminus_+', '+_terminus_-', '-_terminus_+', '-_terminus_-'}
            if u in nodes_to_exclude or v in nodes_to_exclude or G.edges[edge]['branch_point'] in nodes_to_exclude:
                continue

        var_pos = G.get_variant_position(edge)
        is_out_region = len(interval_tree[var_pos]) == 0

        if is_out_region:
            egde_outside_region.append(edge)

    return egde_outside_region

def find_variants_in_intervaltree(interval_tree: IntervalTree,
                                  G: PangenomeGraph,
                                  exclude_terminus: bool = False) -> list:
    egde_in_region = []
    for edge in sorted(list(G.variant_edges)):
        u, v = edge

        if exclude_terminus:
            nodes_to_exclude = {'+_terminus_+', '+_terminus_-', '-_terminus_+', '-_terminus_-'}
            if u in nodes_to_exclude or v in nodes_to_exclude or G.edges[edge]['branch_point'] in nodes_to_exclude:
                continue

        var_pos = G.nodes[u]['position']
        is_in_region = len(interval_tree[var_pos]) > 0

        if is_in_region:
            egde_in_region.append(edge)

    return egde_in_region

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

# Summary functions ----------------------------------------------------------------------------------------------------
def combined_info_vcf_snarl(vcf_path: str, snarl_path: str, level: int = None, pos_range: Tuple[int, int] = None) -> Dict:
    bubble_dict, _ = extract_node_bubble_partition_from_vcf(vcf_path)
    snarl_dict, _ = extract_node_bubble_partition_from_snarl(snarl_path)
    combine_dict = {bubble: bubble_dict.get(bubble, {}) | snarl_dict.get(bubble, {})
                    for bubble in sorted(list(set(bubble_dict.keys()).union(set(snarl_dict.keys()))))}

    if level is not None:
        combine_dict = {key: value for key, value in combine_dict.items() if int(value.get('Level', -1)) == level}
    if pos_range is not None:
        combine_dict = {key: value for key, value in combine_dict.items() if pos_range[0] <= int(value.get('POS', -1)) <= pos_range[1]}

    return combine_dict

def variant_edges_summary(G: PangenomeGraph, var_list: list) -> Dict:
    summary_dict = dict()
    for edge in sorted(list(var_list)):
        if G.is_inversion(edge):
            summary_dict['inversion'] = summary_dict.get('inversion', 0) + 1
        if G.is_replacement(edge):
            summary_dict['replacement'] = summary_dict.get('replacement', 0) + 1
        if G.is_insertion(edge):
            summary_dict['insertion'] = summary_dict.get('insertion', 0) + 1
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

def prepare_dataframe(var_dict):
    return pd.DataFrame({
        "Variant Type": ['Snps', 'Mnps', 'Insertion', 'Deletion', 'Replacement', 'Inversion', 'Back', 'Forward', 'Crossing', 'Total'],
        "Count": [
                  var_dict.get('snps', 0),
                  var_dict.get('mnps', 0),
                  var_dict.get('insertion', 0),
                  var_dict.get('forward_edges', 0),
                  var_dict.get('replacement', 0),
                  var_dict.get('inversion', 0),
                  var_dict.get('back_edges', 0),
                  var_dict.get('forward_edges', 0),
                  var_dict.get('crossing_edges', 0),
                  var_dict.get('total', 0),
        ]
    })

def snarl_level_summary(bubble_dict: Dict) -> Dict:
    level_count = dict()
    for _, bubble_info in bubble_dict.items():
        level_count[int(bubble_info.get('Level', -1))] = level_count.get(int(bubble_info.get('Level', -1)), 0) + 1

    return level_count

def variant_outliers(G: PangenomeGraph,
                     bed_list: list[str],
                     chr_name: str = None,
                     exclude_terminus: bool = False) -> list:
    interval_trees = [get_interval_tree_from_bed(bed, chr_name) for bed in bed_list]
    outside_sets = [set(find_variants_outside_interaltree(interval_tree, G, exclude_terminus)) for interval_tree in
                    interval_trees]
    edge_outside_region = list(set.intersection(*outside_sets))
    return edge_outside_region


def variant_summary_result_by_region(G: PangenomeGraph,
                                    bed_file: str = None,
                                    chr_name: str = None,
                                    exclude_terminus: bool = False) -> Dict:
    if bed_file:
        # def check_edge_in_region(pos_u, pos_v, start, end):
        #     return (start <= pos_u <= end) and (start <= pos_v <= end)
        assert chr_name, "Chromosome name is required for BED file"
        interval_tree = get_interval_tree_from_bed(bed_file, chr_name)
        egde_in_region = find_variants_in_intervaltree(interval_tree, G, exclude_terminus)
    else:
        if exclude_terminus:
            nodes_to_exclude = {'+_terminus_+', '+_terminus_-', '-_terminus_+', '-_terminus_-'}
            egde_in_region = [edge for edge in list(G.variant_edges) if edge[0] not in nodes_to_exclude
                              and edge[1] not in nodes_to_exclude
                              and G.edges[edge]['branch_point'] not in nodes_to_exclude]
        else:
            egde_in_region = G.variant_edges

    var_dict = variant_edges_summary(G, egde_in_region)
    return var_dict

def write_bubble_summary_result(chr_name: str,
                                snarl_path: str,
                                gfa_path: str = None,
                                vcf_path: Optional[str] = None,
                                G_pkl_path: Optional[str] = None,
                                output_dir: Optional[str] = './',
                                method: Optional[str] = "AT"):
    if gfa_path is not None and vcf_path is not None:
        raise ValueError("Both GFA and VCF files are provided. Please provide only one of them.")
    elif gfa_path is None and vcf_path is None:
        raise ValueError("Either GFA or VCF file is required.")

    if not os.path.isfile(snarl_path):
        raise FileNotFoundError(f"Snarl content file not found: {snarl_path}")

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    print("Assigning node to bubbles...")

    if method == "AT":
        bubble_dict, node_partition = extract_node_bubble_partition_from_snarl(snarl_path)
    elif method == "Position":
        raise NotImplementedError("Position method is not implemented yet.")
        # if node_partition is None or bubble_dict is None:
        #     bubble_dict = find_bubbles_from_vcf(G, snarl_path)
        #     node_partition = find_node_partition(G, bubble_dict)
        # method_suffix = "_Position"
    else:
        raise ValueError(f"Invalid method: {method}")

    print("Conducting bubble summary...")
    if gfa_path:
        if not os.path.isfile(gfa_path):
            raise FileNotFoundError(f"GFA file not found: {gfa_path}")

        if G_pkl_path and os.path.isfile(G_pkl_path):
            print(f"Loading graph from {G_pkl_path}...")
            G, walks, walk_sample_names = load_graph_from_pkl(G_pkl_path)
        else:
            print("Constructing graph from GFA...")
            G, walks, walk_sample_names = PangenomeGraph.from_gfa(gfa_path,
                                                                  return_walks=True,
                                                                  compressed=False)

        var_dict_within, var_dict_crossing = get_variants_for_bubbles(G, node_partition)

    if vcf_path:
        if not os.path.isfile(vcf_path):
            raise FileNotFoundError(f"VCF file not found: {vcf_path}")

        var_dict_within, var_dict_crossing = get_variants_for_bubbles_from_vcf(vcf_path, node_partition)


    bubble_var_count_path = f"bubble_variant_counts_{chr_name}_{method}.tsv"

    bubble_list = list(bubble_dict.keys())

    var_with = [var_dict_within.get(key, {}) for key in bubble_list]
    # var_with_summary = [variant_edges_summary(G, var_dict_within.get(key, [])) for key in bubble_list]
    var_cross = [var_dict_crossing.get(key, {}) for key in bubble_list]
    # var_cross_summary = [variant_edges_summary(G, var_dict_crossing.get(key, [])) for key in bubble_list]

    length_with = [len(var_dict_within.get(key, {})) for key in bubble_list]
    length_cross = [len(var_dict_crossing.get(key, {})) for key in bubble_list]
    length_total = [length_with[i] + length_cross[i] for i in range(len(bubble_list))]

    if vcf_path:
        bubble_dict_vcf, _ = extract_node_bubble_partition_from_vcf(vcf_path)
        AC_sum = [sum(ast.literal_eval(f"[{bubble_dict_vcf[x]['AC']}]")) for x in bubble_list]
    else:
        AC_sum = ['.'] * len(bubble_list)

    print("Writing bubble summary to CSV...")
    count_summary_df = pd.DataFrame({
         "Bubble_ids": bubble_list,
         "Total_count": length_total,
         "AC_sum": AC_sum,
         "Within_count": length_with,
         "Crossing_count": length_cross,
         "Within": var_with,
         #"Within_summary": var_with_summary,
         "Crossing": var_cross,
         #"Crossing_summary": var_cross_summary,
         }
    )

    count_summary_df.to_csv(os.path.join(output_dir, bubble_var_count_path), sep='\t')

