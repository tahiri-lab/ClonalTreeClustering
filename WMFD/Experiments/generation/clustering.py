#!/usr/bin/env python3
# -- coding: utf-8 --
"""
clustering.py
-------------
Clustering algorithms for phylogenetic tree distance matrices.
Supports K-Medoids, K-Means, and DBSCAN with evaluation metrics.
"""
import numpy as np
from sklearn.metrics import adjusted_rand_score, silhouette_score, calinski_harabasz_score
from sklearn.manifold import MDS
from sklearn.cluster import KMeans, DBSCAN
# K-Medoids backend
try:
    import kmedoids as KM
    _HAVE_KMEDOIDS = True
except ImportError:
    _HAVE_KMEDOIDS = False

def cluster_data(distance_matrix, method="kmedoids", k=3, **kwargs):
    """
    Main clustering function supporting multiple algorithms.
    
    Args:
        distance_matrix: Precomputed distance matrix (n x n)
        method: "kmedoids", "kmeans", or "dbscan"
        k: Number of clusters (ignored for DBSCAN)
        **kwargs: Additional parameters for specific methods
    Returns:
        labels: Cluster labels array
    """
    D = np.asarray(distance_matrix)
    
    if method == "kmedoids":
        if not _HAVE_KMEDOIDS:
            raise ImportError("Install kmedoids: pip install kmedoids")
        result = KM.fasterpam(D, k)
        return np.array(result.labels, dtype=int)
    
    elif method == "kmeans":
        # Convert distance matrix to coordinates using MDS
        n_components = min(k+1, D.shape[0]-1, 10)  # Reasonable number of dimensions
        mds = MDS(n_components=n_components, dissimilarity="precomputed", 
                  random_state=42, n_init=1, max_iter=300)
        X = mds.fit_transform(D)
        
        kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
        return kmeans.fit_predict(X)
    
    elif method == "dbscan":
        # DBSCAN parameters
        eps = kwargs.get("eps", 0.5)
        min_samples = kwargs.get("min_samples", 5)
        
        dbscan = DBSCAN(eps=eps, min_samples=min_samples, metric="precomputed")
        return dbscan.fit_predict(D)
    
    else:
        raise ValueError(f"Unknown method: {method}. Use 'kmedoids', 'kmeans', or 'dbscan'")

def evaluate_clustering(distance_matrix, labels, y_true=None):
    """
    Evaluate clustering quality with multiple metrics.
    
    Args:
        distance_matrix: Distance matrix used for clustering
        labels: Predicted cluster labels
        y_true: True labels (optional, for ARI calculation)
        
    Returns:
        dict: Evaluation metrics
    """
    results = {}
    D = np.asarray(distance_matrix)
    
    # Silhouette score (works with precomputed distances)
    try:
        if len(set(labels)) > 1:  # Need at least 2 clusters
            results["silhouette"] = silhouette_score(D, labels, metric="precomputed")
        else:
            results["silhouette"] = np.nan
    except Exception:
        results["silhouette"] = np.nan
    
    # Calinski-Harabasz (requires coordinate representation)
    try:
        if len(set(labels)) > 1:
            n_components = min(4, D.shape[0]-1)
            if n_components >= 1:
                mds = MDS(n_components=n_components, dissimilarity="precomputed", 
                         random_state=42, n_init=1, max_iter=200)
                X = mds.fit_transform(D)
                results["calinski_harabasz"] = calinski_harabasz_score(X, labels)
            else:
                results["calinski_harabasz"] = np.nan
        else:
            results["calinski_harabasz"] = np.nan
    except Exception:
        results["calinski_harabasz"] = np.nan
    
    # ARI (if true labels provided)
    if y_true is not None:
        try:
            results["ari"] = adjusted_rand_score(y_true, labels)
        except Exception:
            results["ari"] = np.nan
    
    return results

def process_wmfd_data(wmfd_data, method="kmedoids", **kwargs):
    """
    Apply clustering to WMFD data output.
    
    Args:
        wmfd_data: Output from wmfd.compute_all_wmfd_inmem()
        method: Clustering method to use
        **kwargs: Additional clustering parameters
        
    Returns:
        dict: Results with clustering information
    """
    results = {}
    
    for run_name, entry in wmfd_data.items():
        D = entry["D"]
        true_labels_raw = entry.get("labels_true", [])
        
        # Convert true labels to integers
        if true_labels_raw:
            true_labels = np.array([int(str(lab).split("_")[0]) if isinstance(lab, str) else int(lab) 
                                   for lab in true_labels_raw])
            k = len(np.unique(true_labels))
        else:
            true_labels = None
            k = entry.get("meta", {}).get("K", 3)
        
        # Clustering
        try:
            pred_labels = cluster_data(D, method=method, k=k, **kwargs)
            metrics = evaluate_clustering(D, pred_labels, true_labels)
            
            results[run_name] = {
                "predicted_labels": pred_labels,
                "true_labels": true_labels,
                "metrics": metrics,
                "method": method,
                "k": k
            }
        except Exception as e:
            results[run_name] = {"error": str(e)}
    
    return results

# Test/Demo
if __name__ == "__main__":
    try:
        import wmfd
    except ImportError:
        print("Error: wmfd.py not found")
        exit(1)
    
    print("Generating WMFD data with full parameters...")
    # Utiliser les paramètres par défaut (tous les runs)
    data = wmfd.compute_all_wmfd_inmem(progress=True)
    
    print(f"\nProcessing {len(data)} runs with K-Medoids...")
    results = process_wmfd_data(data, method="kmedoids")
    
    # Afficher tous les résultats
    successful = 0
    ari_scores = []
    
    for run_name, result in results.items():
        if "error" not in result:
            successful += 1
            metrics = result["metrics"]
            ari = metrics.get('ari', 0)
            sil = metrics.get('silhouette', 0)
            ch = metrics.get('calinski_harabasz', 0)
            ari_scores.append(ari)
            
            print(f"{run_name}: ARI={ari:.3f} | Sil={sil:.3f} | CH={ch:.1f}")
        else:
            print(f"{run_name}: ERROR - {result['error']}")
    
    # Résumé statistique
    if ari_scores:
        print(f"\n=== RÉSUMÉ ===")
        print(f"Runs réussis: {successful}/{len(data)}")
        print(f"ARI - Min: {min(ari_scores):.3f}, Max: {max(ari_scores):.3f}, Moyenne: {sum(ari_scores)/len(ari_scores):.3f}")
        
        # Compter les ARI parfaits
        perfect_ari = sum(1 for ari in ari_scores if ari >= 0.999)
        print(f"ARI = 1.000: {perfect_ari}/{len(ari_scores)} runs ({100*perfect_ari/len(ari_scores):.1f}%)")
    
    print(f"\n=== DONNÉES POUR API ===")
    print("Les résultats contiennent:")
    print("- predicted_labels: clusters prédits")  
    print("- true_labels: clusters vrais")
    print("- metrics['ari']: ARI calculé")
    print("Prêt pour l'étape finale API !")