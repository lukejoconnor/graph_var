#!/usr/bin/env python3

from graph import PangenomeGraph

def analyze_gfa(gfa_path, reference_path_index=1):
    """
    Load a GFA file and analyze the edge types.
    
    Args:
        gfa_path (str): Path to the GFA file
        reference_path_index (int): Index of the reference path
    """
    print(f"Loading GFA file: {gfa_path}")
    
    try:
        # Load the graph
        G, walks, walk_sample_names = PangenomeGraph.from_gfa(
            gfa_path,
            return_walks=True,
            compressed=False,
            reference_path_index=reference_path_index
        )
        
        # Get edge type summary
        edge_type_count = G.variant_edges_summary()
        
        print("\nEdge Type Summary:")
        print("-----------------")
        for edge_type, count in edge_type_count.items():
            print(f"{edge_type}: {count}")
            
        print(f"\nTotal number of edges: {G.number_of_edges()}")
        print(f"Total number of nodes: {G.number_of_nodes()}")
        
    except Exception as e:
        print(f"Error analyzing GFA file: {str(e)}")

if __name__ == "__main__":
    # Path to your GFA file
    gfa_path = "data/c4a.gfa"
    analyze_gfa(gfa_path)
