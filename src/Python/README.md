# GBLD PROGRAM

This code comprised of the following main steps:

   - parse_newick: Extracts node weights from Newick strings
   - normalize_weights: Normalizes weights to a 0-1 scale
   - Calculate weight differences: Computes pairwise differences between normalized weights for each pair of trees.
   - calculate_distances: Computes pairwise distances between sequences in a tree
   - normalize_matrix: Normalizes a matrix to a 0-1 scale
   - calculate_W: Calculates the W component of the GBLD metric, then averaging based on the maximum number of nodes
   - calculate_D: Calculates the D component of the GBLD metric, then averaging based on the maximum number of nodes
   - calculate_penalty: Computes the penalty (P) component of the GBLD metric
   - Calculate the final GBLD values.


This code essentially implements the GBLD (Generalized Branch Length Distance) metric for comparing lineage trees, taking into account the overlapping nodes, tree structure, and branch weights.




# newick.from.fasta
For in-depth analysis, we selected a representative file. We assigned unique identifiers to the nucleotide sequences and quantified their occurrence across all files. 

Our analytical process involved several key steps:

1. Sequence Alignment: We utilized the FASTA format and BioPython's AlignIO module to align and read the sequences.

2. Distance Calculation: Employing BioPython's DistanceCalculator, we computed a distance matrix using the 'identity' method.

3. Phylogenetic Tree Construction: We built a phylogenetic tree using the Neighbor-Joining (NJ) method and converted it to Newick format for compatibility with various analysis tools.

4. Frequency and Weighting: We calculated the frequency of each sequence. In the final step, we appended "@" to each sequence name, followed by its abundance value, effectively weighting the sequences.

