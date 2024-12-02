import os
import numpy as np
from itertools import combinations

def read_newick_trees(file_path):
    """
    Read multiple Newick trees from a text file.
    """
    trees = []
    try:
        with open(file_path, 'r') as file:
            for line in file:
                line = line.strip()
                if line:
                    trees.append(line)
    except FileNotFoundError:
        raise FileNotFoundError(f"File not found: {file_path}")
    return trees

def extract_values(newick_str):
    """
    Extract sequence names (including naive), weights, and branch lengths from a Newick string.
    """
    weights = {}
    branch_lengths = {}
    
    parts = newick_str.replace('(', ',').replace(')', ',').strip(';').split(',')
    
    for part in parts:
        part = part.strip()
        if part and '@' in part and ':' in part:
            name_weight, branch = part.split(':')
            name, weight = name_weight.split('@')
            name = name.strip()
            weights[name] = float(weight)
            branch_lengths[name] = float(branch)
    
    return weights, branch_lengths

def sort_sequences(sequences):
    """
    Sort sequences with naive first, then numerically
    """
    def extract_number(seq):
        if seq == 'naive':
            return -1
        return int(seq[3:])
    
    return sorted(sequences, key=extract_number)

def min_max_normalize_matrix(matrix):
    """
    Perform min-max normalization on each column of the matrix.
    """
    normalized = np.zeros_like(matrix, dtype=float)
    
    for j in range(matrix.shape[1]):
        column = matrix[:, j]
        col_min = np.min(column)
        col_max = np.max(column)
        
        if col_max == col_min:
            normalized[:, j] = 0 if col_max == 0 else 1
        else:
            normalized[:, j] = (column - col_min) / (col_max - col_min)
    
    return normalized

def create_matrices(trees):
    """
    Create matrices of weights and branch lengths.
    """
    all_sequences = set()
    for tree in trees:
        weights, _ = extract_values(tree)
        all_sequences.update(weights.keys())
    
    seq_list = sort_sequences(list(all_sequences))
    n_trees = len(trees)
    n_sequences = len(seq_list)
    
    weight_matrix = np.zeros((n_trees, n_sequences))
    branch_matrix = np.zeros((n_trees, n_sequences))
    
    for i, tree in enumerate(trees):
        weights, branches = extract_values(tree)
        for j, seq in enumerate(seq_list):
            weight_matrix[i, j] = weights.get(seq, 0)
            branch_matrix[i, j] = branches.get(seq, 0)
    
    normalized_weight_matrix = min_max_normalize_matrix(weight_matrix)
    normalized_branch_matrix = min_max_normalize_matrix(branch_matrix)
    
    return (normalized_weight_matrix, normalized_branch_matrix, seq_list)

def analyze_tree_pairs(trees):
    """
    Analyze common and uncommon nodes between pairs of trees and calculate normalized sums.
    """
    normalized_weight_matrix, normalized_branch_matrix, seq_list = create_matrices(trees)
    
    weight_results = []
    branch_results = []
    
    total_nodes = len(seq_list)
    
    for (i, tree1), (j, tree2) in combinations(enumerate(trees), 2):
        # Extract nodes for each tree
        nodes1, _ = extract_values(tree1)
        nodes2, _ = extract_values(tree2)
        
        # Find common and uncommon nodes
        common_nodes = set(nodes1.keys()) & set(nodes2.keys())
        all_nodes = set(nodes1.keys()) | set(nodes2.keys())
        uncommon_nodes = all_nodes - common_nodes
        
        num_common = len(common_nodes)
        num_uncommon = len(uncommon_nodes)
        
        # Get indices for common and uncommon nodes
        common_indices = [k for k, seq in enumerate(seq_list) if seq in common_nodes]
        uncommon_indices = [k for k, seq in enumerate(seq_list) if seq in uncommon_nodes]
        
        # Calculate normalized sums for weights
        weight_common_sum = (np.abs(normalized_weight_matrix[i, common_indices] - 
                                  normalized_weight_matrix[j, common_indices]).sum() / 
                           num_common if num_common > 0 else 0)
        
        weight_uncommon_sum = (np.abs(normalized_weight_matrix[i, uncommon_indices] - 
                                    normalized_weight_matrix[j, uncommon_indices]).sum() / 
                             num_uncommon if num_uncommon > 0 else 0)
        
        # Calculate normalized sums for branch lengths
        branch_common_sum = (np.abs(normalized_branch_matrix[i, common_indices] - 
                                  normalized_branch_matrix[j, common_indices]).sum() / 
                           num_common if num_common > 0 else 0)
        
        branch_uncommon_sum = (np.abs(normalized_branch_matrix[i, uncommon_indices] - 
                                    normalized_branch_matrix[j, uncommon_indices]).sum() / 
                             num_uncommon if num_uncommon > 0 else 0)
        
        # Store results with node counts
        weight_results.append({
            'trees': f"Tree_{i+1}_vs_Tree_{j+1}",
            'common_sum': weight_common_sum,
            'uncommon_sum': weight_uncommon_sum,
            'num_common': num_common,
            'num_uncommon': num_uncommon
        })
        
        branch_results.append({
            'trees': f"Tree_{i+1}_vs_Tree_{j+1}",
            'common_sum': branch_common_sum,
            'uncommon_sum': branch_uncommon_sum,
            'num_common': num_common,
            'num_uncommon': num_uncommon
        })
    
    return weight_results, branch_results

def format_results(results, title):
    """
    Format results as a string table with node counts.
    """
    output = [f"\n{title}:", 
              "Tree Pair | Common Nodes Sum | Uncommon Nodes Sum | #Common Nodes | #Uncommon Nodes"]
    output.append("-" * 80)
    
    for result in results:
        output.append(f"{result['trees']:<15} {result['common_sum']:>14.4f} {result['uncommon_sum']:>16.4f} "
                     f"{result['num_common']:>13d} {result['num_uncommon']:>15d}")
    
    return "\n".join(output)

def main():
    base_path = os.path.expanduser("~/1.mahsa.farnia/classificataion_journal")
    input_file = os.path.join(base_path, "weighted_newicks_60.txt")
    output_file = os.path.join(base_path, "node_comparison_normalized_results.txt")
    
    # Read trees and analyze
    trees = read_newick_trees(input_file)
    weight_results, branch_results = analyze_tree_pairs(trees)
    
    # Write results to file
    with open(output_file, 'w') as f:
        f.write(format_results(weight_results, "Normalized Weight Differences (Divided by Node Count)"))
        f.write(format_results(branch_results, "Normalized Branch Length Differences (Divided by Node Count)"))
    
    print(f"Results have been saved to: {output_file}")

if __name__ == "__main__":
    main()