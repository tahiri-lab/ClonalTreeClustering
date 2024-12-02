import os
import pandas as pd
import numpy as np
from sklearn.cluster import DBSCAN

def get_lambda_values():
    """Get lambda values from user input"""
    print("\nPlease enter the lambda values for each metric:")
    lambda1 = float(input("λ₁ (coefficient for Branch Length) = "))
    lambda2 = float(input("λ₂ (coefficient for Weight) = "))
    lambda3 = float(input("λ₃ (coefficient for Degree) = "))
    lambda4 = float(input("λ₄ (coefficient for Height) = "))
    lambda5 = float(input("λ₅ (coefficient for Hamming Distance) = "))
    return lambda1, lambda2, lambda3, lambda4, lambda5

def get_dbscan_params():
    """Get DBSCAN parameters from user input"""
    print("\nPlease enter the DBSCAN parameters:")
    eps = float(input("epsilon (ε) = "))
    min_samples = int(input("minPoints = "))
    return eps, min_samples

def parse_tree_pair(pair_str):
    """Parse tree pair string into two tree names"""
    pair_str = pair_str.strip('()"')
    t1, t2 = pair_str.split(',')
    return t1.strip(), t2.strip()

def calculate_wmfd(row, lambda1, lambda2, lambda3, lambda4, lambda5):
    """Calculate WMFD with debugging information"""
    try:
        penalty = float(row['Penalty'])
        
        common_bl = float(row['Normalized_Common_BL']) if pd.notna(row['Normalized_Common_BL']) else 0
        common_weight = float(row['Normalized_Common_Weight']) if pd.notna(row['Normalized_Common_Weight']) else 0
        common_degree = float(row['Normalized_Common_Degree']) if pd.notna(row['Normalized_Common_Degree']) else 0
        common_height = float(row['Normalized_Common_Height']) if pd.notna(row['Normalized_Common_Height']) else 0
        
        uncommon_bl = float(row['Normalized_Uncommon_BL']) if pd.notna(row['Normalized_Uncommon_BL']) else 0
        uncommon_weight = float(row['Normalized_Uncommon_Weight']) if pd.notna(row['Normalized_Uncommon_Weight']) else 0
        uncommon_degree = float(row['Normalized_Uncommon_Degree']) if pd.notna(row['Normalized_Uncommon_Degree']) else 0
        uncommon_height = float(row['Normalized_Uncommon_Height']) if pd.notna(row['Normalized_Uncommon_Height']) else 0
        
        hamming_dist = float(row['Normalized_Hamming_Distance']) if pd.notna(row['Normalized_Hamming_Distance']) else 0
        
        common_part = (
            lambda1 * common_bl +
            lambda2 * common_weight +
            lambda3 * common_degree +
            lambda4 * common_height
        )
        
        uncommon_part = (
            lambda1 * uncommon_bl +
            lambda2 * uncommon_weight +
            lambda3 * uncommon_degree +
            lambda4 * uncommon_height
        )
        
        wmfd = common_part + (penalty * uncommon_part) + lambda5*hamming_dist
        return wmfd
        
    except Exception as e:
        print(f"Error calculating WMFD for {row['Tree_Pair']}: {str(e)}")
        return None

def create_symmetric_matrix(results_df):
    """Create symmetric matrix from WMFD values"""
    unique_trees = set()
    for pair in results_df['Tree_Pair']:
        t1, t2 = parse_tree_pair(pair)
        unique_trees.add(t1)
        unique_trees.add(t2)
    unique_trees = sorted(list(unique_trees), key=lambda x: int(x.split('_')[1]))
    n = len(unique_trees)
    
    matrix = np.zeros((n, n))
    tree_to_index = {tree: i for i, tree in enumerate(unique_trees)}
    
    for _, row in results_df.iterrows():
        t1, t2 = parse_tree_pair(row['Tree_Pair'])
        i, j = tree_to_index[t1], tree_to_index[t2]
        wmfd_value = row['WMFD']
        if pd.notna(wmfd_value):
            matrix[i, j] = wmfd_value
            matrix[j, i] = wmfd_value
    
    print("\nTree indices:")
    for tree, idx in tree_to_index.items():
        print(f"{tree}: {idx}")
    
    return matrix, unique_trees

def print_symmetric_matrix(matrix, trees):
    """Print symmetric matrix in a readable format"""
    n = len(trees)
    print("\nDistance Matrix:")
    print("-" * 60)
    
    print("Tree", end="")
    for tree in trees:
        print(f"{tree:>8}", end="")
    print()
    
    for i in range(n):
        print(f"{trees[i]:<4}", end="")
        for j in range(n):
            print(f"{matrix[i,j]:8.2f}", end="")
        print()

def main():
    try:
        input_path = os.path.expanduser('~/1.mahsa.farnia/classificataion_journal/tree_metrics 2.csv')
        output_path = os.path.expanduser('~/1.mahsa.farnia/classificataion_journal/wmfd_clustering_results.csv')
        
        print("Reading input file...")
        df = pd.read_csv(input_path)
        
        lambda1, lambda2, lambda3, lambda4, lambda5 = get_lambda_values()
        print("\nCalculating WMFD values...")
        
        results = []
        for idx, row in df.iterrows():
            wmfd = calculate_wmfd(row, lambda1, lambda2, lambda3, lambda4, lambda5)
            results.append({
                'Tree_Pair': row['Tree_Pair'],
                'WMFD': round(wmfd, 4) if wmfd is not None else None
            })
            print(f"WMFD for {row['Tree_Pair']}: {round(wmfd, 4) if wmfd is not None else None}")
        
        results_df = pd.DataFrame(results)
        
        distance_matrix, unique_trees = create_symmetric_matrix(results_df)
        print_symmetric_matrix(distance_matrix, unique_trees)
        
        eps, min_samples = get_dbscan_params()
        
        print("\nDistance Statistics:")
        valid_distances = distance_matrix[distance_matrix > 0]
        print(f"Minimum distance: {np.min(valid_distances):.4f}")
        print(f"Maximum distance: {np.max(valid_distances):.4f}")
        print(f"Mean distance: {np.mean(valid_distances):.4f}")
        
        print("\nPairs within epsilon ({eps}):")
        for i in range(len(unique_trees)):
            for j in range(i+1, len(unique_trees)):
                if distance_matrix[i,j] <= eps and distance_matrix[i,j] > 0:
                    print(f"{unique_trees[i]} - {unique_trees[j]}: {distance_matrix[i,j]:.4f}")
        
        dbscan = DBSCAN(eps=eps, min_samples=min_samples, metric='precomputed')
        labels = dbscan.fit_predict(distance_matrix)
        
        print("\nClustering Results:")
        print("-" * 60)
        n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
        print(f"Number of clusters: {n_clusters}")
        
        for label in sorted(set(labels)):
            if label == -1:
                print("\nNoise points:")
            else:
                print(f"\nCluster {label}:")
            cluster_members = [unique_trees[i] for i, l in enumerate(labels) if l == label]
            print(", ".join(cluster_members))
        
        # Save results
        with open(output_path.replace('.csv', '_results.txt'), 'w') as f:
            f.write("WMFD and Clustering Results\n")
            f.write("=" * 60 + "\n\n")
            
            f.write("Parameters:\n")
            f.write(f"λ₁ (Branch Length) = {lambda1}\n")
            f.write(f"λ₂ (Weight) = {lambda2}\n")
            f.write(f"λ₃ (Degree) = {lambda3}\n")
            f.write(f"λ₄ (Height) = {lambda4}\n")
            f.write(f"λ5 (Height) = {lambda5}\n")
            f.write(f"epsilon (ε) = {eps}\n")
            f.write(f"minPoints = {min_samples}\n\n")
            
            f.write("Distance Matrix:\n")
            n = len(unique_trees)
            f.write("Tree")
            for tree in unique_trees:
                f.write(f"{tree:>8}")
            f.write("\n")
            
            for i in range(n):
                f.write(f"{unique_trees[i]:<4}")
                for j in range(n):
                    f.write(f"{distance_matrix[i,j]:8.2f}")
                f.write("\n")
            
            f.write("\nClustering Results:\n")
            f.write("-" * 40 + "\n")
            for label in sorted(set(labels)):
                if label == -1:
                    f.write("\nNoise points:\n")
                else:
                    f.write(f"\nCluster {label}:\n")
                cluster_members = [unique_trees[i] for i, l in enumerate(labels) if l == label]
                f.write(", ".join(cluster_members) + "\n")
        
        print(f"\nResults have been saved to: {output_path.replace('.csv', '_results.txt')}")
        
    except Exception as e:
        print(f"\nError occurred: {str(e)}")
        raise

if __name__ == "__main__":
    main()