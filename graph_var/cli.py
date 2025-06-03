#!/usr/bin/env python3

import argparse
import sys
import os
from .graph import PangenomeGraph

def parse_args():
    parser = argparse.ArgumentParser(
        description="Identify variants in a pangenome graph and write them to a VCF file",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        "gfa_file",
        help="Input GFA file path"
    )
    
    parser.add_argument(
        "vcf_file",
        help="Output path for VCF file"
    )
    
    parser.add_argument(
        "--chr-id",
        help="Chromosome ID for VCF file (default: chr0)",
        default="chr0"
    )
    
    parser.add_argument(
        "--ref_name",
        help="Name of the reference path (default: 'GRCh38')",
        type=str,
        default='GRCh38'
    )
    
    parser.add_argument(
        "--no-genotypes",
        help="Do not write compute genotypes",
        action="store_true",
        default=False
    )
    
    return parser.parse_args()

def main():
    args = parse_args()
    
    # Check if input file exists
    if not os.path.exists(args.gfa_file):
        print(f"Error: GFA file '{args.gfa_file}' not found", file=sys.stderr)
        sys.exit(1)
    
    # Load the graph
    print(f"Loading GFA file: {args.gfa_file}")
    G = PangenomeGraph.from_gfa_line_by_line(
        args.gfa_file,
        ref_name=args.ref_name
    )
    
    # Generate VCF
    if args.no_genotypes:
        G.write_vcf(None, args.vcf_file, args.chr_id)
    else:
        G.write_vcf(args.gfa_file, args.vcf_file, args.chr_id)

if __name__ == "__main__":
    main()
