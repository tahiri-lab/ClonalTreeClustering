from jaccard_distance_modifie import compute_jaccard_distance_bipartitions
import math

def Threshold(newick_trees, alpha=1, verbose=True, return_threshold=False):
    """
    Computes a frequency threshold for identifying common (frequent) bipartitions across
    a collection of phylogenetic trees, based on the average modified Jaccard distance.

    The threshold is derived by computing the average pairwise Jaccard distance (L)
    between all trees, then scaling it by the number of trees to obtain T = ceil(L * n).
    This value T can be used as a cutoff: a bipartition is considered frequent if it appears 
    in at least T trees.

    Parameters
    ----------
    newick_trees : list of str
        A list of phylogenetic trees represented in Newick format.
    alpha : float, optional (default = 1)
        Weight factor penalizing topological incompatibilities in Jaccard computation.
    verbose : bool, optional (default = True)
        If True, prints the individual distances and summary statistics.
    return_threshold : bool, optional (default = False)
        If True, returns both the average distance (L) and the threshold (T) as a tuple.

    Returns
    -------
    int or tuple(float, int)
        - If return_threshold is False: returns the threshold T (int)
        - If return_threshold is True: returns a tuple (L, T)
          where L is the average distance and T = ceil(L * n)
    """
    n = len(newick_trees)
    if n < 2:
        raise ValueError("At least two trees are required for pairwise comparison.")

    total_distance = 0.0
    pair_count = 0

    # Compute all pairwise modified Jaccard distances
    for i in range(n - 1):
        for j in range(i + 1, n):
            jd = compute_jaccard_distance_bipartitions(newick_trees[i], newick_trees[j], alpha=alpha)
            if verbose:
                print(f"Jaccard Distance (Tree {i+1}, Tree {j+1}) = {jd:.4f}")
            total_distance += jd
            pair_count += 1

    # Compute average distance
    L = total_distance / pair_count

    # Compute frequency threshold (number of trees a bipartition must appear in to be "frequent")
    T = math.ceil(L * n)

    if verbose:
        print(f"\nNumber of trees: {n}")
        print(f"Average Jaccard distance (L): {L:.4f}")
        print(f"Frequency threshold T = ceil(L * n): {T}")

    return (L, T) if return_threshold else T


# ------------------ Example Usage ------------------
if __name__ == "__main__":
    trees = [
        "(((12)2)9,((6,1)7,(3,5,13)10)4,14,(11)15,(16)8)N;",
        "(((2)9)12,((1,4)7,(10,5,13)6)3,14,(15)11,(16)8)N;",
        "(((12)9)2,((13,5)1,(3,4,6)7)10,14,(11)15,(16)8)N;",
        "(((12)2)9,((3,1)7,(13,6,5)4)10,14,(11)15,(8)16)N;",
        "(((9)2)12,((3,1)7,(13,4,6)5)10,14,(11)15,(8)16)N;"
    ]

    threshold = Threshold(trees, alpha=1)
    print(f"\n✅ A bipartition is considered frequent if it appears in ≥ {threshold} trees.")
