# Identify genetic variants in human graph genomes

## Overview
This project, **[Project Title]**, is focused on [Brief Description of the project, including the primary research question or objective]. The main goal of this project is to [Explain the goals or intended outcomes].

## Table of Contents
- [Introduction](#introduction)
- [Installation](#installation)
- [Usage](#usage)

## Introduction
[Provide a detailed introduction to the project, including any relevant background information, key concepts, and why the project is important.]

## Usage
```python
# Example code
from graph import PangenomeGraph

gfa_path = "/path/to/graph.gfa"
G, walks, walk_sample_names = PangenomeGraph.from_gfa(gfa_path, return_walks=True, compressed=False)

# Generate vcf file
vcf_path = "/path/to/output.vcf"
tree_path = "/path/to/output.tree"
chr_id = "chr1"

G.write_vcf(walks, walk_sample_names, vcf_path, tree_path, chr_id, exclude_terminus=True)

```
