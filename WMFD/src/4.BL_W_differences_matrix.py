import os
import numpy as np

def read_newick_trees(file_path):
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
    def extract_number(seq):
        if seq == 'naive':
            return -1
        return int(seq[3:])
    return sorted(sequences, key=extract_number)

def min_max_normalize_matrix(matrix):
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
    
    return normalized_weight_matrix, normalized_branch_matrix, seq_list

def calculate_pairwise_differences(matrix):
    n_trees = matrix.shape[0]
    pairwise_matrices = []
    
    for i in range(n_trees):
        for j in range(i + 1, n_trees):
            differences = np.abs(matrix[i] - matrix[j])
            pairwise_matrices.append((i, j, differences))
    
    return pairwise_matrices

def format_pairwise_differences(pairwise_matrices, sequences, matrix_name):
    lines = []
    lines.append(f"{matrix_name}:")
    header = "Trees      " + "".join(f"{seq:>8}" for seq in sequences)
    lines.append(header)
    lines.append("-" * (10 + 8 * len(sequences)))
    
    for tree1_idx, tree2_idx, differences in pairwise_matrices:
        line = f"T{tree1_idx+1}-T{tree2_idx+1} |"
        line += "".join(f"{val:8.2f}" for val in differences)
        lines.append(line)
    
    lines.append("\n")
    return "\n".join(lines)

def format_matrix_for_file(matrix, sequence_names, matrix_name):
    lines = []
    lines.append(f"{matrix_name}:")
    header = "         " + "".join(f"{seq:>8}" for seq in sequence_names)
    lines.append(header)
    lines.append("-" * (9 + 8 * len(sequence_names)))
    
    for i, row in enumerate(matrix):
        line = f"Tree_{i+1:2} |"
        line += "".join(f"{val:8.2f}" for val in row)
        lines.append(line)
    
    lines.append("\n")
    return "\n".join(lines)

def main():
    base_path = os.path.expanduser("~/1.mahsa.farnia/classificataion_journal")
    input_file = os.path.join(base_path, "weighted_newicks_60.txt")
    output_file = os.path.join(base_path, "normalized_pairwise_differences_60.txt")
    
    trees = read_newick_trees(input_file)
    normalized_weight_matrix, normalized_branch_matrix, sequence_names = create_matrices(trees)
    
    norm_weight_differences = calculate_pairwise_differences(normalized_weight_matrix)
    norm_branch_differences = calculate_pairwise_differences(normalized_branch_matrix)
    
    with open(output_file, 'w') as f:
        f.write("Normalized Matrices:\n")
        f.write(format_matrix_for_file(normalized_weight_matrix, sequence_names, 
                                     "Normalized Weight Matrix"))
        f.write(format_matrix_for_file(normalized_branch_matrix, sequence_names, 
                                     "Normalized Branch Length Matrix"))
        
        f.write("\nNormalized Pairwise Difference Matrices:\n")
        f.write(format_pairwise_differences(norm_weight_differences, sequence_names,
                                          "Normalized Weight Pairwise Differences"))
        f.write(format_pairwise_differences(norm_branch_differences, sequence_names,
                                          "Normalized Branch Length Pairwise Differences"))
    
    print(f"Normalized matrices and their pairwise differences have been saved to: {output_file}")

if __name__ == "__main__":
    main()