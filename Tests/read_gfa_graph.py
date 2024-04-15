from graph import DiED_Graph
from search_tree import max_weight_dfs_tree
import networkx as nx
import pandas as pd

test_graph_1 = DiED_Graph("../grch38_chr20.gfa", reference_path_index=7, spanning_forest_method='dfs')
#test_graph_2 = DiED_Graph("grch38_chr20.gfa", reference_path_index=7, spanning_forest_method='dfs')

snps_dict = test_graph_1.find_snps()

snps_df = pd.DataFrame(
    [(k[0], k[1], v[0], v[1], v[2]) for k, v in snps_dict.items()],
    columns=['SNP_edge_node_v0', 'SNP_edge_node_v1', 'Position', 'REF', 'ALT']
)

snps_df.to_csv("./SNPs.csv")

'''
df1 = pd.read_csv('./chr20_edge_info_0.csv')
df2 = pd.read_csv('./chr20_edge_info_7.csv')
df3 = pd.read_csv('./chr20_edgeinfo.csv')

G = nx.DiGraph()
edges = [
    ('A', 'B', 3), ('A', 'C', 2), ('B', 'D', 4), ('B', 'E', 1),
    ('C', 'F', 5), ('E', 'G', 2), ('F', 'G', 3), ('G', 'H', 1), ('G', 'B', 1),
    ('D', 'H', 2), ('C', 'A', 1), ('F', 'E', 4)
]

G.add_weighted_edges_from(edges)
# Define the reference path (ensure it's valid in the graph)
# reference_path = ['A', 'B', 'E', 'G']
reference_path = ['A', 'B', 'D']
dfs_tree, ref_dag = max_weight_dfs_tree(G, 'A', reference_path)

nx.draw(G)

OutEdgeView([('A', 'B'), ('A', 'C'), ('B', 'E'), ('B', 'D'), ('E', 'G'), ('G', 'H'), ('C', 'F')])
OutEdgeView([('A', 'B'), ('A', 'C'), ('B', 'D'), ('B', 'E'), ('D', 'H'), ('E', 'G'), ('C', 'F')])
'''