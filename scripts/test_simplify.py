from graph_var import PangenomeGraph
import matplotlib.pyplot as plt
import networkx as nx
import os

def main():
    # Read a .gfa file
    gfa_path = "/Users/lukeoconnor/Downloads/grch38_chr20.gfa"
    gfa_path = "data/c4a_with_inversion_and_sequences.gfa"
    G, walks, walk_sample_names = PangenomeGraph.from_gfa(gfa_path,
                                                    return_walks=True, compressed=False)
    
    print(G.variant_edges)
    print(G.ref_alt_alleles(G.variant_edges.pop()))
    simplified_graph = G.delete_small_variants(minimum_allele_length=5)
    print(simplified_graph.number_of_edges(), G.number_of_edges())
    termini = {'+_terminus_+', '-_terminus_+', '+_terminus_-', '-_terminus_-'}
    G.delete_tips(simplified_graph, termini)
    G.contract_paths(simplified_graph, termini)

    # Plotting the graphs
    def plot_graph(G, title, filename):
        plt.figure(figsize=(12, 10))
        pos = nx.spring_layout(G, seed=42)  # For consistent layout between plots
        
        # Draw the graph
        nx.draw(G, pos, with_labels=True, node_color='lightblue', 
                node_size=500, arrows=True, connectionstyle='arc3,rad=0.1')
        
        plt.title(title)
        plt.axis('off')
        plt.tight_layout()
        plt.savefig(filename)
        plt.close()
        print(f"Graph saved to {filename}")
    
    # Create output directory if it doesn't exist
    os.makedirs("output", exist_ok=True)
    
    # Plot original graph
    plot_graph(G, "Original Graph", "output/original_graph.png")
    
    # Plot simplified graph
    plot_graph(simplified_graph, "Simplified Graph", "output/simplified_graph_paths.png")

if __name__ == "__main__":
    main()