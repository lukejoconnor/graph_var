import re

def flip(s):
    if s == '+':
      return '-'
    else:
      return '+'

def read_gfa(filename):
  nodes = []
  sequences = []
  edges = []
  walks = []
  pattern = r'\w+|[<>]'

  with open(filename, 'r') as file:
      for line in file:
          parts = line.strip().split('\t')
          if parts[0] == 'S':
              nodes.append(parts[1])
              sequences.append(parts[2])
          elif parts[0] == 'L':
              edge = (parts[1], parts[3], parts[2], parts[4])
              edges.append(edge)
          elif parts[0] == 'W':
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
  return nodes, edges, walks, sequences
