from graph import PangenomeGraph
import re
import numpy as np
from tqdm import tqdm
from intervaltree import Interval, IntervalTree

def extract_bubble_ids(node_string):
    # Split the string and keep '>' and '<' symbols
    node_ids = re.findall(r'[><]\d+', node_string)

    # Convert the list of strings to a tuple
    node_tuple = tuple(node_ids)
    node_tuple = tuple(map(lambda x: x[1:] + '_+', node_tuple))

    assert len(node_tuple) == 2, f"Invalid bubble ID: {node_string}"

    return node_tuple

def find_bubbles_from_vcf(G: PangenomeGraph, vcf_path: str) -> dict:
    """
    Find bubbles in the pangenome graph that are supported by variants in a VCF file.
    :param pangraph: A PangenomeGraph object
    :param vcf_path: Path to a VCF file
    :return: A list of tuples representing the bubbles found
    """
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

def find_vertice_partition(G: PangenomeGraph, bubble_dict: dict) -> dict:
    """
    Find the partition of the pangenome graph based on the bubbles found.
    :param pangraph: A PangenomeGraph object
    :param bubble_dict: A dictionary of bubbles found
    :return: A list of tuples representing the partition
    """
    node_partition = {}
    position_sorted_nodes = sorted(G.nodes(data=True), key=lambda x: x[1]['position'])
    # Create the interval tree
    position_interval_tree = IntervalTree(Interval(min(v[0], v[1]), max(v[0], v[1]), k) for k, v in bubble_dict.items())

    for node, node_attr in position_sorted_nodes:
        node_partition[node] = [interval.data for interval in position_interval_tree[node_attr['position']]]

    return node_partition


def count_num_variants_in_bubbles(G: PangenomeGraph, bubble_dict: dict) -> dict:
    """
    Count the number of variants in each bubble.
    :param pangraph: A PangenomeGraph object
    :param bubble_dict: A dictionary of bubbles found
    :param vcf_path: Path to a VCF file
    :return: A dictionary containing the number of variants in each bubble
    """
    bubble_variants = {}

    for k, v in bubble_dict.items():
        start, end = min(v[0], v[1]), max(v[0], v[1])
        bubble_variants[k] = len(G.get_variants_at_interval((start, end)))

    return bubble_variants





