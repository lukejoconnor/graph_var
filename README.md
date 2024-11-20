# Identifying genetic variants in the pangenome using a reference tree

## Table of Contents
- [Introduction](#introduction)
- [Installation](#installation)
- [Usage](#usage)
- [Command Line Interface](#command-line-interface)

## Introduction
This repository is very much a work in progress. If you use our code, please do reach out with questions and feedback.

## Usage
```python
from graph import PangenomeGraph

# Read a .gfa file
gfa_path = "/path/to/graph.gfa"
reference_path_index = 1
G, walks, walk_sample_names = PangenomeGraph.from_gfa(gfa_path, 
                                                return_walks=True, compressed=False, 
                                                reference_path_index=reference_path_index)

# Generate vcf file
vcf_path = "/path/to/output.vcf"
tree_path = "/path/to/output.tree"
chr_id = "chr1"
G.write_vcf(walks, walk_sample_names, vcf_path, tree_path, chr_id, exclude_terminus=True)

# Enumerate variants of different types
edge_type_count: dict = G.variant_edges_summary()

# Get the genotype of a walk, then reconstruct edge visit counts
genotype: dict = G.genotype(walks[0])
edge_visit_counts: dict = G.count_edge_visits(genotype)

```

## Command Line Interface

```bash
graph_var <gfa_file> [options]
```

### Required Arguments
- `gfa_file`: Path to the input GFA file containing the pangenome graph

### Optional Arguments
- `--nodeinfo`: Path to input node info file (CSV format)
- `--edgeinfo`: Path to input edge info file (CSV format)
- `--out-nodeinfo`: Output path for node info file (CSV format)
- `--out-edgeinfo`: Output path for edge info file (CSV format)
- `--vcf`: Output path for VCF file
- `--chr-id`: Chromosome ID for VCF file (default: "chr0")
- `--ref-path-index`: Index of the reference path (default: None, will try to find GRCh38)
- `--summary [file]`: Print summary of variant types. If a file path is provided, write the summary to that file instead of stdout

### Example Usage
```bash
# Basic usage - analyze a GFA file and output VCF
# Will automatically try to find GRCh38 reference path
graph_var input.gfa --vcf output.vcf

# Generate node and edge info files (CSV format)
graph_var input.gfa --out-nodeinfo nodes.csv --out-edgeinfo edges.csv

# Analyze with specific reference path and chromosome ID
graph_var input.gfa --vcf output.vcf --ref-path-index 2 --chr-id chr20

# Print variant summary to console
graph_var input.gfa --summary

# Write variant summary to a file
graph_var input.gfa --summary summary.txt
