from graph import PangenomeGraph

def main():
    # Read a .gfa file
    gfa_path = "data/c4a_with_inversion.gfa"
    reference_path_index = 1
    G, walks, walk_sample_names = PangenomeGraph.from_gfa(gfa_path,
                                                    return_walks=True, compressed=False)
    # Generate vcf file
    vcf_path = "c4a_with_inversion.vcf"
    tree_path = "c4a_with_inversion.tree"
    chr_id = "chr6"
    G.write_vcf(walks, walk_sample_names, vcf_path, tree_path, chr_id, exclude_terminus=True)

    # Enumerate variants of different types
    edge_type_count: dict = G.variant_edges_summary()
    print(edge_type_count)

    # Get the genotype of a walk, then reconstruct edge visit counts
    genotype: dict = G.genotype(walks[5])
    edge_visit_counts: dict = G.count_edge_visits(genotype)
    # print(walks[5], edge_visit_counts)

if __name__ == "__main__":
    main()