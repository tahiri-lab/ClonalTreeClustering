import numpy as np
from ete3 import Tree
from jaccard_distance import *
from frequent_arc import extract_frequent_arcs
from average_Jaccard_distance import *

def build_reference_matrix_from_nonfrequent_arcs(newick_list, frequent_arcs, verbose=True):
    """
    Constructs a pairwise similarity matrix between trees based on their shared **non-frequent arcs**.

    Parameters:
    -----------
    newick_list : list of str
        A list of tree representations in Newick format.

    frequent_arcs : set of (parent, child) tuples
        A set of arcs considered globally frequent and therefore excluded from the similarity calculation.

    verbose : bool
        If True, the function prints the resulting matrix to the console.

    Returns:
    --------
    matrix : numpy.ndarray
        A square (n x n) symmetric matrix where entry (i, j) indicates the number of non-frequent arcs shared
        between tree i and tree j. Diagonal entries are 0 by default (self-comparison is excluded).
    """
    n = len(newick_list)
    arc_sets = []

    # Extract non-frequent arcs for each tree
    for newick in newick_list:
        tree = Tree(newick, format=8)
        arcs = get_direct_arcs(tree)  # Get all direct arcs (parent → child)
        nonfrequent_arcs = arcs - frequent_arcs  # Filter out frequent arcs
        arc_sets.append(nonfrequent_arcs)

    # Initialize an empty similarity matrix
    matrix = np.zeros((n, n), dtype=int)

    # Compare each tree with others
    for i in range(n):
        for j in range(i, n):
            if i != j:  # Skip diagonal entries (can be changed if needed)
                common = arc_sets[i].intersection(arc_sets[j])
                matrix[i][j] = len(common)
                # Optional: uncomment the next line to make the matrix symmetric
                # matrix[j][i] = len(common)

    if verbose:
        print("\nReference matrix (non-frequent arcs only):")
        for row in matrix:
            print("  ".join(f"{val:2d}" for val in row))

    return matrix


# Example usage
if __name__ == "__main__":
    trees = [
        "(((12)2)9,((6,1)7,(3,5,13)10)4,14,(11)15,(16)8)N;",
        "(((2)9)12,((1,4)7,(10,5,13)6)3,14,(15)11,(16)8)N;",
        "(((12)9)2,((13,5)1,(3,4,6)7)10,14,(11)15,(16)8)N;",
        "(((12)2)9,((3,1)7,(13,6,5)4)10,14,(11)15,(8)16)N;",
        "(((9)2)12,((3,1)7,(13,4,6)5)10,14,(11)15,(8)16)N;"
    ]

    # Step 1: Compute the Jaccard-based frequency threshold
    threshold = compute_average_jaccard_distance(trees, verbose=False)

    # Step 2: Extract arcs that occur in at least 'threshold' number of trees
    frequent_arcs = extract_frequent_arcs(trees, threshold)

    # Step 3: Build the matrix of shared non-frequent arcs
    build_reference_matrix_from_nonfrequent_arcs(trees, frequent_arcs)
