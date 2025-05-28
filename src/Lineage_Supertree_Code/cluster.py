from Reference_tab import *

def cluster_finder(reference_tab):
    """
    Identifies clusters of similar elements based on their maximum similarity values
    from a reference similarity matrix.

    Parameters:
    reference_tab (np.ndarray): A square 2D numpy array representing similarity scores 
                                between pairs of trees.

    Returns:
    list of list of int: A list where each element is a list of indices representing a cluster. 
                         Each cluster groups indices whose similarity with the reference row 
                         reaches the maximum similarity value for that row.
    """
    all_clusters = []
    all_max = np.max(reference_tab, axis=1)  # Get the max similarity for each row
    n = len(reference_tab)

    for i in range(n - 1):
        cluster = [i]
        max_sim = all_max[i]
        for j in range(i + 1, n):
            # Add index j to the cluster if its similarity with i equals i's max similarity
            if reference_tab[i][j] == max_sim:
                cluster.append(j)
        all_clusters.append(cluster)

    return all_clusters


# Example usage
if __name__ == "__main__":
    # A set of Newick-formatted phylogenetic trees
    trees = [
        "(((12)2)9,((6,1)7,(3,5,13)10)4,14,(11)15,(16)8)N;",
        "(((2)9)12,((1,4)7,(10,5,13)6)3,14,(15)11,(16)8)N;",
        "(((12)9)2,((13,5)1,(3,4,6)7)10,14,(11)15,(16)8)N;",
        "(((12)2)9,((3,1)7,(13,6,5)4)10,14,(11)15,(8)16)N;",
        "(((9)2)12,((3,1)7,(13,4,6)5)10,14,(11)15,(8)16)N;"
    ]
    
    # Step 1: Compute the average Jaccard distance to determine the threshold
    threshold = compute_average_jaccard_distance(trees, verbose=False)

    # Step 2: Extract arcs (tree structures) that occur in at least 'threshold' number of trees
    frequent_arcs = extract_frequent_arcs(trees, threshold, verbose=False)

    # Step 3: Construct a similarity matrix using arcs *not* identified as frequent
    ref_matrix = build_reference_matrix_from_nonfrequent_arcs(trees, frequent_arcs, verbose=False)

    # Step 4: Group trees into clusters based on their highest similarity values
    clusters = cluster_finder(ref_matrix)

    # Output the list of clusters
    print(clusters)
