import numpy as np
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, calinski_harabasz_score
from Bipartition_matrix import build_bipartition_matrix

def analyse_kmeans_bipartitions(matrix, k_range=None):
    """
    Performs K-Means clustering analysis on a binary matrix representing 
    presence/absence of bipartitions across trees.

    Parameters
    ----------
    matrix : list or 2D np.array
        Binary matrix (trees × bipartitions), where each row corresponds to a tree 
        and each column to a bipartition. The last row (bipartition frequencies) 
        is excluded from clustering.

    k_range : iterable of int, optional
        Range of cluster numbers to test (e.g., range(2, n_trees)). 
        If None, defaults to 2 to n_trees - 1.

    Displays
    -------
    - Clustering metrics (Silhouette score, Calinski-Harabasz index, Inertia)
    - Cluster assignments for each tree
    """
    # Remove the last row (frequency row)
    matrix_sans_freq = matrix[:-1]
    arbres = [f"Tree {i+1}" for i in range(len(matrix_sans_freq))]

    X = np.array(matrix_sans_freq)
    n_trees = X.shape[0]

    # Default range of clusters: from 2 to n_trees - 1
    if k_range is None:
        k_range = range(2, n_trees)

    for k in k_range:
        print(f"\n===== k = {k} clusters (K-Means) =====")
        kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
        labels = kmeans.fit_predict(X)

        # Compute clustering quality metrics
        sil = silhouette_score(X, labels)
        ch = calinski_harabasz_score(X, labels)
        inertia = kmeans.inertia_

        print(f"Silhouette Score        : {sil:.2f}")
        print(f"Calinski-Harabasz Index : {ch:.2f}")
        print(f"Inertia (within-cluster sum of squares): {inertia:.2f}")

        # Show members of each cluster
        for cluster_id in range(k):
            members = [arbres[i] for i in range(n_trees) if labels[i] == cluster_id]
            print(f" - Cluster {cluster_id + 1} : {', '.join(members)}")
     

# ------------------------- Example Usage -------------------------
if __name__ == "__main__":

    trees_newick = [
        "(((12)2)9,((1,6)7,(3,5,13)10)4,14,(11)15,(16)8)N;",
        "(((2)9)12,((1,4)7,(10,5,13)6)3,14,(15)11,(16)8)N;",
        "(((12)9)2,((13,5)1,(3,4,6)7)10,14,(11)15,(16)8)N;",
        "(((12)2)9,((1,3)7,(5,6,13)4)10,14,(11)15,(8)16)N;",
        "(((9)2)12,((3,1)7,(13,4,6)5)10,14,(11)15,(8)16)N;"
    ]

    # Construct bipartition matrix from input trees
    matrix, bipartitions = build_bipartition_matrix(trees_newick)

    # Perform K-Means clustering on tree-bipartition profiles
    analyse_kmeans_bipartitions(matrix)
