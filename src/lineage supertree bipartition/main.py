from bipart_of_one_tree import *
from Bipartition_matrix import *
from newick_reader import *
from common_structure import *
from number_incompatibility import *
from threshold import *
from cluster_finder_Kmean import *
from cluster_finder_Cmean import *

# -----------------------------
# Example Newick-formatted trees
# -----------------------------
trees_newick = [
    "(((12)2)9,((1,6)7,(3,5,13)10)4,14,(11)15,(16)8)N;",
    "(((2)9)12,((1,4)7,(10,5,13)6)3,14,(15)11,(16)8)N;",
    "(((12)9)2,((13,5)1,(3,4,6)7)10,14,(11)15,(16)8)N;",
    "(((12)2)9,((1,3)7,(5,6,13)4)10,14,(11)15,(8)16)N;",
    "(((9)2)12,((3,1)7,(13,4,6)5)10,14,(11)15,(8)16)N;"
]

# -----------------------------
# Build the bipartition matrix
# -----------------------------
matrix, bipartitions = build_bipartition_matrix(trees_newick)

# Display the binary matrix (trees × bipartitions)
print("Binary Tree × Bipartition Matrix:\n")
for i, row in enumerate(matrix, 1):
    if i == len(trees_newick) + 1:
        print(f"Frequency: {row}")
    else:
        print(f"Tree {i}: {row}")

# Display bipartitions with indices
print("\nBipartition Columns:")
for i, bipart in enumerate(bipartitions, 1):
    print(f"{i:02d}: {sorted(bipart)}")

# -----------------------------
# Compute threshold value
# -----------------------------
threshold_value = Threshold(trees_newick, alpha=1, verbose=False)
print(f"\n✅ A bipartition is considered frequent if it appears in ≥ {threshold_value} trees.")

# -----------------------------
# Extract common structure
# -----------------------------
common_structure_of_all_trees = Common_structure(trees_newick, bipartitions, matrix, threshold_value)
print("✅ Frequent Bipartitions (Common Structure):")
for b in common_structure_of_all_trees:
    print(b)

# -----------------------------
# Run clustering analyses
# -----------------------------
analyse_kmeans_bipartitions(matrix)
analyse_cmeans_bipartitions(matrix)
