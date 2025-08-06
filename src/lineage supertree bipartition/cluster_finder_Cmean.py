import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import silhouette_score
from fcmeans import FCM
from sklearn.decomposition import PCA
from Bipartition_matrix import build_bipartition_matrix

def analyse_cmeans_bipartitions(matrix, k_range=None):
    """
    Performs fuzzy clustering (Fuzzy C-Means) on a binary tree-bipartition matrix.

    Parameters
    ----------
    matrix : list or 2D np.array
        Binary matrix (trees × bipartitions), with or without the last frequency row.
    
    k_range : iterable of int, optional
        Range of cluster numbers to test (e.g., range(2, len(trees))).
        Defaults to 2 to (n_trees - 1).

    Displays
    -------
    - Silhouette score for each cluster number (based on hard assignments)
    - Degree of membership for each tree to each cluster
    - Heatmap of membership degrees
    - PCA projection for visual clustering assessment
    """

    # Remove frequency row (last row)
    matrix_sans_freq = matrix[:-1]
    tree_labels = [f"Tree {i+1}" for i in range(len(matrix_sans_freq))]

    X = np.array(matrix_sans_freq)
    n_trees = X.shape[0]

    # Default cluster range: from 2 to n_trees - 1
    if k_range is None:
        k_range = range(2, n_trees)

    for k in k_range:
        print(f"\n===== k = {k} clusters (Fuzzy C-Means) =====")
        fcm = FCM(n_clusters=k, random_state=42)
        fcm.fit(X)

        u_matrix = fcm.u  # Membership matrix
        hard_labels = np.argmax(u_matrix, axis=1)  # Hard cluster assignment for silhouette

        sil_score = silhouette_score(X, hard_labels)
        print(f"Silhouette Score (based on hard labels): {sil_score:.2f}")

        # Display membership degrees for each tree
        for i, tree in enumerate(tree_labels):
            memberships = [f"C{j+1}:{u_matrix[i][j]:.2f}" for j in range(k)]
            print(f" - {tree} : {', '.join(memberships)}")

        # Heatmap of membership degrees
        import pandas as pd
        membership_df = pd.DataFrame(u_matrix, columns=[f"C{i+1}" for i in range(k)], index=tree_labels)
        plt.figure(figsize=(8, 4))
        sns.heatmap(membership_df, annot=True, cmap="YlGnBu", cbar=True)
        plt.title(f"Fuzzy C-Means - Membership Degrees (k={k})")
        plt.tight_layout()
        plt.show()

    # PCA projection for visual representation
    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X)

    for k in k_range:
        fcm = FCM(n_clusters=k, random_state=42)
        fcm.fit(X)
        hard_labels = np.argmax(fcm.u, axis=1)

        plt.figure(figsize=(6, 5))
        for cluster_id in range(k):
            indices = np.where(hard_labels == cluster_id)[0]
            plt.scatter(X_pca[indices, 0], X_pca[indices, 1], label=f"Cluster {cluster_id + 1}")
        
        # Annotate trees
        for i, txt in enumerate(tree_labels):
            plt.annotate(txt, (X_pca[i, 0], X_pca[i, 1]))

        plt.title(f"Fuzzy C-Means (k={k}) - PCA Projection")
        plt.legend()
        plt.tight_layout()
        plt.show()


# ------------------------- Example Usage -------------------------
if __name__ == "__main__":

    trees_newick = [
        "(((12)2)9,((1,6)7,(3,5,13)10)4,14,(11)15,(16)8)N;",
        "(((2)9)12,((1,4)7,(10,5,13)6)3,14,(15)11,(16)8)N;",
        "(((12)9)2,((13,5)1,(3,4,6)7)10,14,(11)15,(16)8)N;",
        "(((12)2)9,((1,3)7,(5,6,13)4)10,14,(11)15,(8)16)N;",
        "(((9)2)12,((3,1)7,(13,4,6)5)10,14,(11)15,(8)16)N;"
    ]

    matrix, bipartitions = build_bipartition_matrix(trees_newick)

    analyse_cmeans_bipartitions(matrix)
