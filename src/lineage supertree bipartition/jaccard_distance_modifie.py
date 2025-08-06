from bipart_of_one_tree import extract_bipartitions
from number_incompatibility import number_of_incompatibility

def compute_jaccard_distance_bipartitions(newick1, newick2, alpha=1):
    """
    Computes a modified Jaccard distance between two phylogenetic trees in Newick format,
    incorporating a penalty based on topological incompatibility between bipartitions.

    This method extends the classical Jaccard index by subtracting a scaled incompatibility
    score from the intersection cardinality. This allows for a finer-grained measure of 
    dissimilarity that accounts for conflicting splits.

    Parameters
    ----------
    newick1 : str
        The first phylogenetic tree in Newick format.
    newick2 : str
        The second phylogenetic tree in Newick format.
    alpha : float, optional (default = 1)
        The penalty weight applied to each incompatibility between bipartitions.

    Returns
    -------
    float
        A normalized Jaccard distance between 0 and 1:
        - 0 indicates identical trees (no distance),
        - 1 indicates complete dissimilarity.
    """

    # Extract sets of bipartitions from both trees
    bipart_1 = {frozenset(bp) for bp in extract_bipartitions(newick1)}
    bipart_2 = {frozenset(bp) for bp in extract_bipartitions(newick2)}

    # Compute set intersection and union
    intersection = bipart_1.intersection(bipart_2)
    union = bipart_1.union(bipart_2)

    # Compute the number of incompatible bipartitions
    incompatibility = number_of_incompatibility(newick1, newick2)

    # Handle edge case: empty bipartitions (e.g., empty or trivial trees)
    if len(union) == 0:
        return 0.0

    # Modified Jaccard distance: penalized overlap
    penalized_similarity = max(0, len(intersection) - (alpha * incompatibility))
    jaccard_distance = 1 - penalized_similarity / len(union)

    return jaccard_distance


# ------------------ Example Usage ------------------
if __name__ == "__main__":
    # Example trees in Newick format
    newick_str1 = "(((12)2)9,((1,6)7,(3,5,13)10)4,14,(11)15,(16)8)N;"
    newick_str2 = "(((2)9)12,((1,4)7,(10,5,13)6)3,14,(15)11,(16)8)N;"
    newick_str3 = "(((12)9)2,((13,5)1,(3,4,6)7)10,14,(11)15,(16)8)N;"
    newick_str4 = "(((12)2)9,((1,3)7,(5,6,13)4)10,14,(11)15,(8)16)N;"
    newick_str5 = "(((9)2)12,((3,1)7,(13,4,6)5)10,14,(11)15,(8)16)N;"

    # Compute Jaccard distance with default alpha
    distance_default = compute_jaccard_distance_bipartitions(newick_str2, newick_str3)
    print(f"Jaccard distance (default alpha=1): {distance_default:.4f}")

    # Compute Jaccard distance with custom penalty weight
    distance_custom = compute_jaccard_distance_bipartitions(newick_str2, newick_str3, alpha=0.5)
    print(f"Jaccard distance (custom alpha=0.5): {distance_custom:.4f}")
