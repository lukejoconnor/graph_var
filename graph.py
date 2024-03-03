from math import inf
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from utils import flip, read_gfa
from search_tree import max_weight_dfs_tree

class DiED_Graph(nx.DiGraph):
    def __init__(self, gfa_file=None):
        super(DiED_Graph, self).__init__()
        self.reference_tree = nx.DiGraph()
        self.path_basis = None
        self.walks = None
        self.reference_dag = nx.DiGraph()
        self.reference_path = None
        self.position = [{} for i in range(2)]

        if gfa_file:
          genes, edges, walks = read_gfa(gfa_file)
          print("Num of Genes:", len(genes))
          print("Num of Edges:", len(edges))
          for gene in genes:
            self.add_gene(gene)
          for edge in edges:
            self.add_diedge(edge[0], edge[1], edge[2], edge[3])
            self.add_diedge(edge[1], edge[0], flip(edge[3]), flip(edge[2]))
          self.walks = walks

        # assign weights to edges
        n_missing = 0
        n_nonmissing = 0
        for walk in self.walks:
          for i in range(len(walk) - 1):
              u, v = walk[i], walk[i + 1]
              if not self.has_edge(u, v):  # Check if the edge exists in the graph
                n_missing += 1
                if n_missing <= 10:
                  print(u,v)
                continue
              #self[u][v]['weight'] += 1  # Increase the weight of the edge by 1
              n_nonmissing += 1
        print(n_missing, n_nonmissing)

        # create a DFS forest
        source_nodes = [node for node, in_degree in self.in_degree() if in_degree == 0]
        sink_nodes = [node for node, out_degree in self.out_degree() if out_degree == 0]

        self.start_node = 'start_node'  # This is one way to ensure the new node has a unique ID
        self.add_node(self.start_node)
        self.end_node = 'end_node'  # This is one way to ensure the new node has a unique ID
        self.add_node(self.end_node)

        for source_node in source_nodes:
            self.add_edge(self.start_node, source_node, weight=0)
        for sink_node in sink_nodes:
            self.add_edge(sink_node, self.end_node, weight=0)

        self.add_weights()

        # self.linear_reference = self.walks[0]

        # Construct reference tree and ref DAG
        self.reference_tree, self.reference_dag = max_weight_dfs_tree(self, self.start_node) #TODO: maybe replace this function

        # Ref path: from start to end node in the ref dag
        self.reference_path = nx.shortest_path(self.reference_tree, self.start_node, self.end_node)

        # Define node positions
        for direction in range(2):
          print(direction)
          self.get_node_positions(direction)

        # Check consistency of node positions
        for node in self.reference_dag.nodes():
          assert(self.position[0][node] <= self.position[1][node])

        # All of the edges in the pangenome graph not in the reference tree
        all_edges = set(self.edges())
        reference_edges = set(self.reference_tree.edges())
        self.path_basis = all_edges.difference(reference_edges)


    def add_gene(self, gene, seq=''):

      self.add_node(str(gene)+'_+_i', sequence=seq, direction='+')
      self.add_node(str(gene)+'_+_o', sequence=seq, direction='+')
      self.add_node(str(gene)+'_-_i', sequence=seq, direction='-')
      self.add_node(str(gene)+'_-_o', sequence=seq, direction='-')
      self.add_edge(str(gene)+'_+_i', str(gene)+'_+_o', kind='intra', weight=0)
      self.add_edge(str(gene)+'_-_i', str(gene)+'_-_o', kind='intra', weight=0)

    def add_diedge(self, gene1, gene2, direction1, direction2):

      self.add_edge(str(gene1)+'_'+direction1+'_o', str(gene2)+'_'+direction2+'_i', kind='inter', weight=0)

    def add_weights(self):
      for walk_list in self.walks:
        for i in range(int(len(walk_list))-1):
          self.edges[walk_list[i], walk_list[i+1]]['weight'] += 1

    def add_weights_test(self):
      for id, walk_list in enumerate(self.walks):
        for i in range(int(len(walk_list))-1):
          if id == 0:
            self.edges[walk_list[i], walk_list[i+1]]['weight'] += 1e+9
          else:
            self.edges[walk_list[i], walk_list[i+1]]['weight'] += 1

    def get_node_positions(self, direction):
      print(direction)
      self.position[direction] = {u:-inf * (-1)**direction for u in self.reference_dag.nodes}
      for n, u in enumerate(self.reference_path):
        self.position[direction][u] = n
      assert(self.position[direction][self.start_node] == 0)
      assert(self.position[direction][self.end_node] > 0)

      if direction == 0:
        G = self.reference_dag
      else:
        G = self.reference_dag.reverse()

      order = nx.topological_sort(G)
      for u in order:
        pred = list(G.predecessors(u))
        pred.append(u)

        if direction == 0:
          # print(u, self.position[direction][u])
          self.position[direction][u] = np.max([self.position[direction][v] for v in pred])
          # print(u, self.position[direction][u])
        else:
          self.position[direction][u] = np.min([self.position[direction][v] for v in pred])


    def find_path_to_linear_ref(self, bfs_tree, node):
      current = node
      while current not in self.linear_reference + [self.start_node]:
          # In a BFS tree, each node except the root has exactly one predecessor
          preds = list(bfs_tree.predecessors(current))
          current = preds[0]
      return current

    def show(self, G, i_pos=None):
      if i_pos == None:
        pos = nx.spring_layout(G)
      else:
        pos = i_pos

      plt.figure(figsize=(15,10))
      # Draw the graph
      nx.draw_networkx_nodes(G, pos, node_color='white', node_shape='o', edgecolors='black', node_size=500)

      # Draw edges
      for edge, key in enumerate(G.edges):
        if key[2] == 'r':
          nx.draw_networkx_edges(G, pos, arrows=True, arrowstyle='-', edgelist=[key], connectionstyle=f'arc3,rad={0}', edge_color='red')
        else:
          nx.draw_networkx_edges(G, pos, arrows=True, arrowstyle='-', edgelist=[key], connectionstyle=f'arc3,rad={0.30}', edge_color='black')

      # Draw node labels
      nx.draw_networkx_labels(G, pos, font_size=8)

      # Show plot
      plt.show()
