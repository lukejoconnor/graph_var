from math import inf
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from utils import flip, read_gfa
from search_tree import max_weight_dfs_tree
from kruskal import kruskal_mst
from icecream import ic

class DiED_Graph(nx.DiGraph):
    def __init__(self,
                 gfa_file: str = None,
                 reference_walk_index: int = None,
                 spanning_forest_method: str = 'dfs'
                 ):
        super(DiED_Graph, self).__init__()
        self.reference_tree: nx.classes.digraph.DiGraph = nx.DiGraph()
        self.reference_dag: nx.classes.digraph.DiGraph = nx.DiGraph()
        self.variant_edges: set = set()
        self.walks: list = []
        self.reference_path: list = []
        self.position: dict = {}
        self.num_walks: int = 0
        self.num_variants: int = 0
        self.genotypes: list = []

        if gfa_file:
            genes, edges, walks = read_gfa(gfa_file)
            print("Num of Genes:", len(genes))
            print("Num of Edges:", len(edges))
            for gene in genes:
                self.add_gene(gene)
            for n, edge in enumerate(edges):
                self.add_diedge(edge[0], edge[1], edge[2], edge[3], n)
                self.add_diedge(edge[1], edge[0], flip(edge[3]), flip(edge[2]), -n)
            self.walks = walks
            self.num_walks = len(walks)

        # Add universal source and sink nodes
        self.add_source_and_sink_nodes()

        # Add weights to edges
        self.add_weights()

        # Create spanning tree
        self.add_reference(reference_walk_index= reference_walk_index, spanning_forest_method = spanning_forest_method)

        # Call variants
        self.count_variants()

    def add_source_and_sink_nodes(self):
        source_nodes = [node for node, in_degree in self.in_degree() if in_degree == 0]
        sink_nodes = [node for node, out_degree in self.out_degree() if out_degree == 0]

        self.start_node = 'start_node'
        self.add_node(self.start_node)
        self.end_node = 'end_node'
        self.add_node(self.end_node)

        for source_node in source_nodes:
            self.add_edge(self.start_node, source_node, weight=0)
        for sink_node in sink_nodes:
            self.add_edge(sink_node, self.end_node, weight=0)

    def add_reference(self, reference_walk_index: int = None, spanning_forest_method: str = 'dfs'):

        initial_reference_tree = nx.DiGraph()
        if reference_walk_index is not None:
            if reference_walk_index >= self.num_walks or reference_walk_index < 0:
                raise ValueError(f'Reference walk index should be an integer >= 0 and < {self.num_walks}')
            self.reference_path = self.walks[reference_walk_index]

        # Construct reference tree and ref DAG
        if spanning_forest_method == 'kruskal':
            for u, v in zip(self.reference_path, self.reference_path[1:]):
                # Step 3: Check if the edge exists in G and then add it to H with the same attributes
                if self.has_edge(u, v):
                    edge_attr = self[u][v]  # Get the attributes of the edge from G
                    initial_reference_tree.add_edge(u, v, **edge_attr)  # Add the edge to H with the attributes

            self.reference_tree = kruskal_mst(
                self,
                universal_source=self.start_node,
                initial_graph=initial_reference_tree
            )
        elif spanning_forest_method == 'dfs':
            self.reference_tree, self.reference_dag = max_weight_dfs_tree(self,
                                                                          source=self.start_node,
                                                                          reference_path=self.reference_path)
        else:
            raise ValueError(f"spanning_forest_method should be either 'kruskal' or 'dfs'")


        # All of the edges in the pangenome graph not in the reference tree
        all_edges = set(self.edges())
        reference_edges = set(self.reference_tree.edges())
        self.variant_edges = all_edges.difference(reference_edges)
        self.num_variants = len(self.variant_edges)

    def add_positions(self):
        if self.reference_dag.number_of_nodes() < self.number_of_nodes():
            raise ValueError('Reference DAG does not have enough nodes, probably because it is not defined')

        # Define node positions
        for direction in range(2):
            self.get_node_positions(direction)

        # Check consistency of node positions
        for node in self.reference_dag.nodes():
            assert (self.position[0][node] <= self.position[1][node])

    def add_gene(self, gene, seq=''):

      self.add_node(str(gene)+'_+', sequence=seq, direction='+')
      self.add_node(str(gene)+'_-', sequence=seq, direction='-')

    def add_diedge(self, gene1, gene2, direction1, direction2, n):

      self.add_edge(str(gene1)+'_'+direction1, str(gene2)+'_'+direction2, weight=0, index=n)

    def add_weights(self):
      for walk_list in self.walks:
        for i in range(int(len(walk_list))-1):
          self.edges[walk_list[i], walk_list[i+1]]['weight'] += 1

    def count_variants(self):

        self.genotypes = [{} for i in range(self.num_walks)]
        for individual in range(self.num_walks):
            walk_list = self.walks[individual]
            for i in range(len(walk_list) - 1):
                edge = walk_list[i], walk_list[i + 1]
                if edge not in self.variant_edges:
                    continue

                if edge in self.genotypes[individual]:
                    self.genotypes[individual][edge] += 1
                else:
                    self.genotypes[individual][edge] = 1

    def get_node_positions(self, direction):
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
