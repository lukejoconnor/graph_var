from graph import PangenomeGraph
import re
import numpy as np
from tqdm import tqdm
from intervaltree import Interval, IntervalTree
from typing import List, Tuple, Union, Dict, Set, Optional
from collections import defaultdict
from utils import _node_convert

def write_dfs_tree_to_gfa(G: PangenomeGraph, filename: str):
    with open(filename, 'wb') as gfa_file:
        for node, node_attr in G.reference_tree.nodes(data=True):
            node_id, symbol = node.split('_')
            gfa_file.write(f'S\t{node_id}\t{node_attr["sequence"]}\n'.encode())

        for edge in G.reference_tree.edges(data=True):
            u, v, edge_attr = edge
            u_id, u_symbol = u.split('_')
            v_id, v_symbol = v.split('_')
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
            AT = [x for x in INFO_list if x.startswith('AT=')][0]
            bubble_id_tuple = extract_bubble_ids(ID)
            nodes_list = extract_nodes_in_bubble(AT)
            bubbles[bubble_id_tuple] = nodes_list
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
