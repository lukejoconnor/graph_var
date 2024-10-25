# Identifying genetic variants in the pangenome using a reference tree

## Table of Contents
- [Introduction](#introduction)
- [Installation](#installation)
- [Usage](#usage)

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
