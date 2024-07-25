# GBLD PROGRAM

This code comprised of the following main steps:

   - parse_newick: Extracts node weights from Newick strings
   - normalize_weights: Normalizes weights to a 0-1 scale
   - Calculate differences between normalized weights for each pair of trees.
   - calculate_distances: Computes pairwise distances between sequences in a tree
   - normalize_matrix: Normalizes a matrix to a 0-1 scale
   - calculate_W: Calculates the W component of the GBLD metric, then averaging based on the maximum number of nodes
   - calculate_D: Calculates the D component of the GBLD metric, then averaging based on the maximum number of nodes
   - calculate_penalty: Computes the penalty (P) component of the GBLD metric
   - Calculate the final GBLD values.


This code essentially implements the GBLD (Generalized Branch Length Distance) metric for comparing lineage trees, taking into account the overlapping nodes, tree structure, and branch weights.
