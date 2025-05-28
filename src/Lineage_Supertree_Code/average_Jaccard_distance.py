from jaccard_distance import *
import math

def compute_average_jaccard_distance(newick_trees, verbose=True, return_threshold=False):
    """
    Computes the average pairwise Jaccard distance (L) across a list of phylogenetic trees
    provided in Newick format, and determines a threshold frequency value for identifying
    frequent arcs based on the formula T = ceil(L * n), where n is the number of trees.

    Parameters:
    ----------
    newick_trees : list of str
        A list of phylogenetic trees in Newick format.

    verbose : bool, optional (default=True)
        If True, prints summary information including the number of trees, average Jaccard
        distance, and the frequency threshold.

    return_threshold : bool, optional (default=False)
        If True, returns both the average Jaccard distance (L) and the computed threshold T.
        Otherwise, returns only the threshold T.

    Returns:
    -------
    int or tuple
        - If return_threshold=False: returns T = ceil(L * n), the minimal frequency threshold
          for defining common arcs.
        - If return_threshold=True: returns a tuple (L, T).
    """
    n = len(newick_trees)
    if n < 2:
        raise ValueError("At least two trees are required for pairwise comparison.")

    total_distance = 0.0
    for i in range(n - 1):
        for j in range(i + 1, n):
            jd = compute_jaccard_distance_arcs(newick_trees[i], newick_trees[j])
            total_distance += jd

    L = (2 / (n * (n - 1))) * total_distance
    T = math.ceil(L * n)

    if verbose:
        print(f"Number of trees: {n}")
        print(f"Average Jaccard distance (L): {L:.4f}")
        print(f"Threshold (ceil(L * n)): {T}")

    return (L, T) if return_threshold else T


# Example usage :
if __name__ == "__main__":
    trees = [
        "(((12)2)9,((6,1)7,(3,5,13)10)4,14,(11)15,(16)8)N;",
        "(((2)9)12,((1,4)7,(10,5,13)6)3,14,(15)11,(16)8)N;",
        "(((12)9)2,((13,5)1,(3,4,6)7)10,14,(11)15,(16)8)N;",
        "(((12)2)9,((3,1)7,(13,6,5)4)10,14,(11)15,(8)16)N;",
        "(((9)2)12,((3,1)7,(13,4,6)5)10,14,(11)15,(8)16)N;"
    ]

    threshold = compute_average_jaccard_distance(trees)
    print(f"An arc is considered frequent if it appears in ≥ {threshold} trees.")
