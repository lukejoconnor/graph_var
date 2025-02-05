import re
import gzip
import pickle

def sequence_complement(s: str) -> str:
    def _base_complement(letter: str) -> str:
        base_dict: dict = {'A': 'T', 'C': 'G', 'T' : 'A', 'G': 'C'}
        return base_dict[letter] if letter in base_dict else letter

    return ''.join(_base_complement(c) for c in reversed(s))

def node_complement(s: str):
    return s[:-1] + _flip(s[-1])

def edge_complement(e: tuple[str, str]):
    return node_complement(e[1]), node_complement(e[0])

def walk_complement(w: list[str]) -> list[str]:
    return [node_complement(v) for v in w[::-1]]

def _flip(s):
    if s == '+':
        return '-'
    elif s == '-':
        return '+'
    else:
        raise ValueError()

def _node_recover(node_id):
    node, direction = node_id.split('_')
    if direction == '+':
        return '>'+node
    elif direction == '-':
        return '<'+node
    else:
        raise ValueError()

def _node_convert(node_id):
    direction, node = node_id[0], node_id[1:]
    if direction == '>':
        return node+"_+"
    elif direction == '<':
        return node+"_-"
    else:
        raise ValueError()

def read_gfa(filename, compressed=False):
  nodes = []
  sequences = []
  edges = []
  walks = []
  walk_sample_names = []
  pattern = r'\w+|[<>]'
  reference_index = 0
  hit_reference = False

  data_dict = {}

  if compressed:
      file = gzip.open(filename, 'rt')
  else:
      file = open(filename, 'r')
  for line in file:
      parts = line.strip().split('\t')
      if parts[0] == 'S':
          nodes.append(parts[1])
          sequences.append(parts[2])
      elif parts[0] == 'L':
          edge = (parts[1], parts[3], parts[2], parts[4])
          edges.append(edge)
      elif parts[0] == 'W':
          if not hit_reference:
              if parts[1] == 'GRCh38':
                  hit_reference = True
              else:
                  reference_index += 1
          sample_name = parts[1]+'_'+parts[2]
          p = parts[6]
          matches = re.findall(pattern, p)
          # List to store node IDs
          node_ids = []

          # Iterate through the matches to generate node IDs
          for i, match in enumerate(matches):
              if match not in ['<', '>']:  # If the match is a word
                  if matches[i - 1] == '<':
                      # For '<', generate IDs with '-'
                      node_ids.append(f'{match}_-')
                  elif matches[i - 1] == '>':
                      # For '>', generate IDs with '+'
                      node_ids.append(f'{match}_+')
          #node_ids.append('end_node')
          walks.append(node_ids)
          walk_sample_names.append(sample_name)

  data_dict['nodes'] = nodes
  data_dict['edges'] = edges
  data_dict['walks'] = walks
  data_dict['walk_sample_names'] = walk_sample_names
  data_dict['sequences'] = sequences
  data_dict['reference_index'] = reference_index

  return data_dict

def save_graph_to_pkl(G, walks, walk_sample_names, path, compressed=False):
    if compressed:
        with gzip.open(path, 'wb') as file:
            pickle.dump(G, file, protocol=pickle.HIGHEST_PROTOCOL)
            pickle.dump(walks, file, protocol=pickle.HIGHEST_PROTOCOL)
            pickle.dump(walk_sample_names, file, protocol=pickle.HIGHEST_PROTOCOL)
    else:
        with open(path, 'wb') as file:
            pickle.dump(G, file, protocol=pickle.HIGHEST_PROTOCOL)
            pickle.dump(walks, file, protocol=pickle.HIGHEST_PROTOCOL)
            pickle.dump(walk_sample_names, file, protocol=pickle.HIGHEST_PROTOCOL)

def load_graph_from_pkl(path, compressed=False):
    if compressed:
        with gzip.open(path, 'rb') as file:
            G = pickle.load(file)
            walks = pickle.load(file)
            walk_sample_names = pickle.load(file)
    else:
        with open(path, 'rb') as file:
            G = pickle.load(file)
            walks = pickle.load(file)
            walk_sample_names = pickle.load(file)
    return G, walks, walk_sample_names

def group_walks_by_name(walks: list, names: list) -> dict:
    """
    Group walks by their names into sublists.
    
    Args:
        walks: List of walks (each walk is a list of nodes)
        names: List of names corresponding to each walk
    
    Returns:
        List of lists, where each sublist contains all walks with the same name
    """
    # Create a dictionary to store walks by name
    from collections import defaultdict
    walks_by_name = defaultdict(list)
    
    # Group walks by their names
    for walk, name in zip(walks, names):
        walks_by_name[name].append(walk)
    
    return walks_by_name

def nearly_identical_alleles(allele1: str, allele2: str):
    """
    Returns True if two sequences differ by at most one base.

    Args:
        allele1: String representation of an allele
        allele2: String representation of an allele

    Returns:
        True if the alleles are nearly identical, False otherwise
    """
    length_difference = len(allele1) - len(allele2)
    if abs(length_difference) > 1:
        return False
    mismatches = 0
    idx1, idx2 = 0, 0
    while idx1 < len(allele1) and idx2 < len(allele2):
        base1 = allele1[idx1]
        base2 = allele2[idx2]
        if base1 == base2:
            idx1 += 1
            idx2 += 1
            continue
        mismatches += 1
        if mismatches > 1:
            return False
        if length_difference >= 0:  # SNP or deletion
            idx1 += 1
        if length_difference <= 0:  # SNP or insertion
            idx2 += 1

    return True