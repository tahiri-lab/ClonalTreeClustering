from bipart_of_one_tree import *

def build_bipartition_matrix(list_of_newicks):
    """
    Constructs a binary bipartition matrix from a list of phylogenetic trees.

    Each row (except the last) corresponds to a tree and each column to a unique bipartition.
    A value of 1 indicates the presence of a bipartition in a tree, 0 indicates its absence.
    The final row contains the frequency of each bipartition across all input trees.

    Parameters
    ----------
    list_of_newicks : List[str]
        List of trees in Newick format.

    Returns
    -------
    matrix : List[List[int]]
        Binary matrix representing bipartition presence/absence per tree, 
        followed by a row of bipartition frequencies.
        
    all_bipartitions : List[frozenset]
        Ordered list of all unique bipartitions found across all trees.
    """

    # A set to collect all unique bipartitions found across all trees
    all_bipartitions = set()

    # A list to store the set of bipartitions for each individual tree
    bipartitions_per_tree = []

    # Extract bipartitions from each tree and update the global set
    for newick in list_of_newicks:
        biparts = extract_bipartitions(newick)
        biparts_set = set(biparts)
        bipartitions_per_tree.append(biparts_set)
        all_bipartitions.update(biparts_set)

    # Sort bipartitions for consistent column ordering across runs
    all_bipartitions = sorted(all_bipartitions, key=lambda x: (len(x), sorted(x)))

    # Construct the binary presence/absence matrix
    matrix = []
    for biparts_set in bipartitions_per_tree:
        row = [1 if bipart in biparts_set else 0 for bipart in all_bipartitions]
        matrix.append(row)

    # Compute and append a frequency row: number of trees each bipartition appears in
    frequencies = [sum(row[i] for row in matrix) for i in range(len(all_bipartitions))]
    matrix.append(frequencies)

    return matrix, all_bipartitions


# ----------------- Example Usage -----------------
if __name__ == "__main__":
    # Input: a list of trees in Newick format
    trees_newick = [
        "(((12)2)9,((1,6)7,(3,5,13)10)4,14,(11)15,(16)8)N;",
        "(((2)9)12,((1,4)7,(10,5,13)6)3,14,(15)11,(16)8)N;",
        "(((12)9)2,((13,5)1,(3,4,6)7)10,14,(11)15,(16)8)N;",
        "(((12)2)9,((1,3)7,(5,6,13)4)10,14,(11)15,(8)16)N;",
        "(((9)2)12,((3,1)7,(13,4,6)5)10,14,(11)15,(8)16)N;"
    ]

    # Build the bipartition matrix
    matrix, bipartitions = build_bipartition_matrix(trees_newick)

    # Print the binary matrix (one row per tree, last row is frequency)
    for i, row in enumerate(matrix, 1):
        print(f"Tree {i}: {row}")

    # Display the list of all bipartitions with their index and taxa
    print("\nBipartition columns:")
    for i, bipart in enumerate(bipartitions, 1):
        print(f"{i:02d}: {sorted(bipart)}")
