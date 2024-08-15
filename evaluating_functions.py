from graph import PangenomeGraph
import re
import numpy as np
from tqdm import tqdm
from intervaltree import Interval, IntervalTree
from typing import List, Tuple, Union, Dict, Set, Optional
from collections import defaultdict


def extract_bubble_ids(node_string: str) -> Tuple:
    # Split the string and keep '>' and '<' symbols
    node_ids = re.findall(r'[><]\d+', node_string)

    # Convert the list of strings to a tuple
    node_tuple = tuple(map(lambda x: x[1:] + '_+', node_ids))

    assert len(node_tuple) == 2, f"Invalid bubble ID: {node_string}"

    return node_tuple


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


def find_node_partition(G: PangenomeGraph, bubble_dict: Dict) -> Tuple[Dict, Dict]:
    # Small epsilon value to shift the start of intervals
    epsilon = 1e-1

    node_partition_u, node_partition_v = {}, {}
    position_sorted_nodes = sorted(G.nodes(data=True), key=lambda x: x[1]['position'])
    # Create the interval tree
    left_close_position_interval_tree = IntervalTree(Interval(min(v), max(v), k) for k, v in bubble_dict.items())
    right_close_position_interval_tree = IntervalTree(Interval(min(v) + epsilon, max(v) + epsilon, k) for k, v in bubble_dict.items())

    for node, node_attr in position_sorted_nodes:
        node_partition_u[node] = set([interval.data for interval in left_close_position_interval_tree[node_attr['position']]])
        node_partition_v[node] = set([interval.data for interval in right_close_position_interval_tree[node_attr['position']]])

    return (node_partition_u, node_partition_v)


def get_variants_for_bubbles(G: PangenomeGraph,
                             node_partition_u: Dict[str, Set[Tuple[str, str]]],
                             node_partition_v: Dict[str, Set[Tuple[str, str]]],
                             inversion_only: bool = False,
                             ) -> Tuple[Dict, Dict, Dict]:
    bubble_within_variants = defaultdict(set)
    bubble_crossing_variants = defaultdict(set)
    bubble_all_variants = defaultdict(set)

    for edge in G.variant_edges:
        u, v = edge

        if inversion_only:
            if not G.is_inversion(edge):
                continue

        bubbles_u = node_partition_u[u]
        bubbles_v = node_partition_v[v]

        bubble_set_union = bubbles_u.union(bubbles_v)
        bubble_set_intersection = bubbles_u.intersection(bubbles_v)
        bubble_set_complementary = bubble_set_union - bubble_set_intersection

        for bubble in bubble_set_intersection:
            bubble_within_variants[bubble].add(edge)

        for bubble in bubble_set_complementary:
            bubble_crossing_variants[bubble].add(edge)

        for bubble in bubble_set_union:
            bubble_all_variants[bubble].add(edge)

    return bubble_within_variants, bubble_crossing_variants, bubble_all_variants


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
