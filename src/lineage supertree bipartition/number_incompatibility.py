from bipart_of_one_tree import extract_bipartitions

def number_of_incompatibility(tree1, tree2):
    """
    Computes the number of incompatibilities between the bipartitions of two phylogenetic trees.

    Two bipartitions A and B are considered incompatible if they overlap but neither is a subset of the other.
    This is based on the classical definition of bipartition incompatibility in consensus and supertree analysis.

    Parameters
    ----------
    tree1 : str
        A phylogenetic tree in Newick format.
    tree2 : str
        Another phylogenetic tree in Newick format.

    Returns
    -------
    int
        The number of pairwise incompatible bipartitions between the two trees.
    """

    nb_incompatibility = 0

    # Extract bipartitions from each input tree
    bipartitions_1 = extract_bipartitions(tree1)
    bipartitions_2 = extract_bipartitions(tree2)

    # Compare each bipartition from tree1 with each from tree2
    for A in bipartitions_1:
        A_set = set(A)
        for B in bipartitions_2:
            B_set = set(B)

            # A and B are incompatible if:
            # - They are not disjoint
            # - A is not a subset of B
            # - B is not a subset of A
            if not (A_set & B_set == set() or A_set <= B_set or B_set <= A_set):
                nb_incompatibility += 1

    return nb_incompatibility


# ----------------- Example Usage -----------------
if __name__ == "__main__":
    # Define several trees in Newick format
    newick_str1 = "(((12)2)9,((1,6)7,(3,5,13)10)4,14,(11)15,(16)8)N;"
    newick_str2 = "(((2)9)12,((1,4)7,(10,5,13)6)3,14,(15)11,(16)8)N;"
    newick_str3 = "(((12)9)2,((13,5)1,(3,4,6)7)10,14,(11)15,(16)8)N;"
    newick_str4 = "(((12)2)9,((1,3)7,(5,6,13)4)10,14,(11)15,(8)16)N;"
    newick_str5 = "(((9)2)12,((3,1)7,(13,4,6)5)10,14,(11)15,(8)16)N;"

    # Compute the incompatibility count between two trees
    result = number_of_incompatibility(newick_str1, newick_str2)
    print("Number of incompatibilities:", result)
