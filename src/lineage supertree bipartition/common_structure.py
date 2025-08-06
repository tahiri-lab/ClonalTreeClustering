from threshold import Threshold
from Bipartition_matrix import build_bipartition_matrix

def Common_structure(trees, bipartitions, matrix, threshold_value):
    """
    Identifies the common structure across a set of phylogenetic trees by 
    extracting frequent bipartitions that occur in at least 'threshold_value' trees.

    Parameters
    ----------
    trees : list of str
        List of phylogenetic trees in Newick format.

    bipartitions : list of frozenset
        The full set of unique bipartitions observed across all trees.

    matrix : list of list of int
        Binary presence/absence matrix of shape (n_trees + 1, n_bipartitions), 
        where each row corresponds to a tree and each column to a bipartition.
        The last row contains the frequency (number of trees) in which 
        each bipartition appears.

    threshold_value : int
        Minimum number of trees in which a bipartition must appear to be 
        considered frequent (i.e., part of the common structure).

    Returns
    -------
    list of frozenset
        A list of bipartitions that represent the common structure, 
        i.e., bipartitions that appear in at least `threshold_value` trees.
    """

    bipart_of_common_structure = []

    # The last row of the matrix contains the frequency of each bipartition
    freq_row = matrix[-1]

    # Iterate over all bipartitions and retain those that are frequent
    for i, freq in enumerate(freq_row):
        if freq >= threshold_value:
            bipart_of_common_structure.append(bipartitions[i])

    return bipart_of_common_structure


# --------------------- Example Usage ---------------------
if __name__ == "__main__":
    trees = [
        "(((12)2)9,((6,1)7,(3,5,13)10)4,14,(11)15,(16)8)N;",
        "(((2)9)12,((1,4)7,(10,5,13)6)3,14,(15)11,(16)8)N;",
        "(((12)9)2,((13,5)1,(3,4,6)7)10,14,(11)15,(16)8)N;",
        "(((12)2)9,((3,1)7,(13,6,5)4)10,14,(11)15,(8)16)N;",
        "(((9)2)12,((3,1)7,(13,4,6)5)10,14,(11)15,(8)16)N;"
    ]
    
    # Compute the frequency threshold based on average modified Jaccard distances
    threshold_value = Threshold(trees)
    
    # Build the bipartition presence/absence matrix
    matrix, biparts = build_bipartition_matrix(trees)
    
    # Extract the bipartitions common to at least 'threshold_value' trees
    common_biparts = Common_structure(trees, biparts, matrix, threshold_value)

    # Display frequent bipartitions representing the consensus structure
    print("✅ Frequent bipartitions (common structure):")
    for b in common_biparts:
        print(b)
