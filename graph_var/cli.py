#!/usr/bin/env python3

import argparse
import sys
import os
import tempfile
from .graph import PangenomeGraph
import json

def parse_args():
    parser = argparse.ArgumentParser(
        description="Analyze genetic variants in pangenome graphs",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        "gfa_file",
        help="Input GFA file path"
    )
    
    parser.add_argument(
        "--nodeinfo",
        help="Input node info file (optional)"
    )
    
    parser.add_argument(
        "--edgeinfo",
        help="Input edge info file (optional)"
    )
    
    parser.add_argument(
        "--out-nodeinfo",
        help="Output path for node info file"
    )
    
    parser.add_argument(
        "--out-edgeinfo",
        help="Output path for edge info file"
    )
    
    parser.add_argument(
        "--vcf",
        help="Output path for VCF file"
    )
    
    parser.add_argument(
        "--chr-id",
        help="Chromosome ID for VCF file (default: chr0)",
        default="chr0"
    )
    
    parser.add_argument(
        "--ref-path-index",
        help="Index of the reference path (default: None, will try to find GRCh38)",
        type=int,
        default=None
    )
    
    parser.add_argument(
        "--summary",
        nargs='?',
        const=True,
        help="Print summary of variant types. Optionally specify a file to write the summary to"
    )
    
    return parser.parse_args()

def write_info_file(data, filepath):
    """Write node or edge info to a JSON file"""
    with open(filepath, 'w') as f:
        json.dump(data, f, indent=2)

def main():
    args = parse_args()
    
    # Check if input file exists
    if not os.path.exists(args.gfa_file):
        print(f"Error: GFA file '{args.gfa_file}' not found", file=sys.stderr)
        sys.exit(1)
    
    # Create output directories if they don't exist
    output_files = [args.vcf, args.out_nodeinfo, args.out_edgeinfo]
    for filepath in output_files:
        if filepath:
            os.makedirs(os.path.dirname(os.path.abspath(filepath)), exist_ok=True)
    
    # Load the graph
    try:
        print(f"Loading GFA file: {args.gfa_file}")
        G, walks, walk_sample_names = PangenomeGraph.from_gfa(
            args.gfa_file,
            return_walks=True,
            compressed=False,
            reference_path_index=args.ref_path_index
        )
        
        # Write node info if requested
        if args.out_nodeinfo:
            G.write_nodeinfo(args.out_nodeinfo)
            print(f"Node info written to: {args.out_nodeinfo}")
        
        # Write edge info if requested
        if args.out_edgeinfo:
            G.write_edgeinfo(args.out_edgeinfo)
            print(f"Edge info written to: {args.out_edgeinfo}")
        
        # Generate VCF if requested
        if args.vcf:
            # Create a temporary file for the tree
            with tempfile.NamedTemporaryFile(suffix='.tree', delete=False) as tmp_tree:
                tree_path = tmp_tree.name
                try:
                    G.write_vcf(walks, walk_sample_names, args.vcf, tree_path, args.chr_id)
                    print(f"VCF file written to: {args.vcf}")
                finally:
                    # Clean up the temporary tree file
                    if os.path.exists(tree_path):
                        os.unlink(tree_path)
        
        # Print summary if requested
        if args.summary:
            summary = G.variant_edges_summary()
            summary_text = "\nVariant Types Summary:\n" + "\n".join(f"{k}: {v}" for k, v in summary.items())
            
            if isinstance(args.summary, str):
                # Write to file if a filename was provided
                os.makedirs(os.path.dirname(os.path.abspath(args.summary)), exist_ok=True)
                with open(args.summary, 'w') as f:
                    f.write(summary_text)
                print(f"Summary written to: {args.summary}")
            else:
                # Print to stdout if no filename was provided
                print(summary_text)
            
            print(f"\nGraph Statistics:")
            print(f"Total nodes: {G.number_of_nodes()}")
            print(f"Total edges: {G.number_of_edges()}")
            print(f"Number of walks: {len(walks)}")
        
    except Exception as e:
        print(f"Error processing GFA file: {str(e)}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
