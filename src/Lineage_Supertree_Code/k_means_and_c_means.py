import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, calinski_harabasz_score, davies_bouldin_score
from fcmeans import FCM
import numpy as np

# === 1. Chargement des données ===
fichier_csv = "tableau bipart.csv"  # Chemin vers ton fichier
df = pd.read_csv(fichier_csv, sep=';', index_col=0)
X = df.values
arbres = df.index.tolist()
k_range = range(2, min(len(df), 5))

# === 2. Clustering K-Means (dur) ===
for k in k_range:
    print(f"\n===== k = {k} clusters (K-Means) =====")
    kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
    labels = kmeans.fit_predict(X)
    
    sil = silhouette_score(X, labels)
    ch = calinski_harabasz_score(X, labels)
    db = davies_bouldin_score(X, labels)
    inertia = kmeans.inertia_

    print(f"Silhouette : {sil:.2f}")
    print(f"Calinski-Harabasz : {ch:.2f}")
    print(f"Davies-Bouldin : {db:.2f}")
    print(f"Inertie : {inertia:.2f}")

    for cluster_id in range(k):
        membres = [arbres[i] for i in range(len(arbres)) if labels[i] == cluster_id]
        print(f" - Cluster {cluster_id + 1} : {', '.join(membres)}")

# === 3. Clustering flou (Fuzzy C-Means) ===
for k in k_range:
    print(f"\n===== k = {k} clusters (Fuzzy C-Means) =====")
    fcm = FCM(n_clusters=k, random_state=42)
    fcm.fit(X)
    u_matrix = fcm.u
    hard_labels = np.argmax(u_matrix, axis=1)
    
    sil_fuzzy = silhouette_score(X, hard_labels)
    print(f"Silhouette (dur) : {sil_fuzzy:.2f}")

    # Détails des degrés d'appartenance
    for i, arbre in enumerate(arbres):
        parts = [f"C{j+1}:{u_matrix[i][j]:.2f}" for j in range(k)]
        print(f" - {arbre} : {', '.join(parts)}")

    # === 4. Heatmap ===
    membership_df = pd.DataFrame(u_matrix, columns=[f"C{i+1}" for i in range(k)], index=arbres)
    plt.figure(figsize=(8, 4))
    sns.heatmap(membership_df, annot=True, cmap="YlGnBu", cbar=True)
    plt.title(f"Fuzzy C-Means - Degrés d'appartenance (k={k})")
    plt.tight_layout()
    plt.show()

# === 5. PCA projection (2D) ===
pca = PCA(n_components=2)
X_pca = pca.fit_transform(X)

# K-Means projections
for k in k_range:
    kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
    labels = kmeans.fit_predict(X)

    plt.figure(figsize=(6, 5))
    for cluster_id in range(k):
        indices = np.where(labels == cluster_id)
        plt.scatter(X_pca[indices, 0], X_pca[indices, 1], label=f"Cluster {cluster_id+1}")
    for i, txt in enumerate(arbres):
        plt.annotate(txt, (X_pca[i, 0], X_pca[i, 1]))
    plt.title(f"K-Means (k={k}) - PCA projection")
    plt.legend()
    plt.tight_layout()
    plt.show()

# FCM projections
for k in k_range:
    fcm = FCM(n_clusters=k, random_state=42)
    fcm.fit(X)
    hard_labels = np.argmax(fcm.u, axis=1)

    plt.figure(figsize=(6, 5))
    for cluster_id in range(k):
        indices = np.where(hard_labels == cluster_id)
        plt.scatter(X_pca[indices, 0], X_pca[indices, 1], label=f"Cluster {cluster_id+1}")
    for i, txt in enumerate(arbres):
        plt.annotate(txt, (X_pca[i, 0], X_pca[i, 1]))
    plt.title(f"Fuzzy C-Means (k={k}) - PCA projection")
    plt.legend()
    plt.tight_layout()
    plt.show()

