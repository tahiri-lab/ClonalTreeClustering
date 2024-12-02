import os
import re
import numpy as np

def get_node_heights(newick_str):
    if ': ' in newick_str:
        newick_str = newick_str.split(': ')[1]
    
    heights = {'naive': 1}
    
    def get_next_node(s, start):
        i = start
        parentheses = 0
        while i < len(s):
            if s[i] == '(':
                parentheses += 1
            elif s[i] == ')':
                parentheses -= 1
            elif s[i] == ',' and parentheses == 0:
                break
            i += 1
        return s[start:i], i
    
    current_height = 2
    to_process = [(newick_str[1:-2], 'naive')]
    
    while to_process:
        next_level = []
        for content, parent in to_process:
            i = 0
            while i < len(content):
                if content[i] == '(':
                    count = 1
                    j = i + 1
                    while count > 0:
                        if content[j] == '(': count += 1
                        if content[j] == ')': count -= 1
                        j += 1
                    
                    k = j
                    while k < len(content) and content[k] not in ',);':
                        k += 1
                    
                    if j < k:
                        node = content[j:k].split('@')[0].strip()
                        heights[node] = current_height
                        next_level.append((content[i+1:j-1], node))
                    i = k
                else:
                    node_str, next_i = get_next_node(content, i)
                    if node_str:
                        node = node_str.split('@')[0].strip()
                        heights[node] = current_height
                    i = next_i + 1
        
        to_process = next_level
        current_height += 1
    
    return heights

def create_height_matrix(trees):
    all_nodes = set()
    tree_heights = []
    
    for tree in trees:
        heights = get_node_heights(tree)
        tree_heights.append(heights)
        all_nodes.update(heights.keys())
    
    sorted_nodes = sorted(all_nodes, key=lambda x: (0 if x == 'naive' else 
                                                  1 if x.startswith('seq') else 2,
                                                  int(re.search(r'\d+', x).group()) if re.search(r'\d+', x) else 0))
    
    height_matrix = []
    for heights in tree_heights:
        row = [heights.get(node, 0) for node in sorted_nodes]
        height_matrix.append(row)
    
    return np.array(height_matrix), sorted_nodes

def normalize_matrix(matrix):
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

def calculate_differences(normalized_matrix):
    n_trees = normalized_matrix.shape[0]
    differences = []
    
    for i in range(n_trees):
        for j in range(i + 1, n_trees):
            diff = np.abs(normalized_matrix[i] - normalized_matrix[j])
            differences.append((f"Tree_{i+1}-Tree_{j+1}", diff))
    
    return differences

def print_matrices(height_matrix, normalized_matrix, differences, nodes, output_file):
    with open(output_file, 'w') as f:
        # Original Height Matrix
        f.write("Original Height Matrix:\n\n")
        header = "         " + "".join(f"{node:>8}" for node in nodes)
        f.write(header + "\n")
        f.write("-" * (9 + 8 * len(nodes)) + "\n")
        
        for i, row in enumerate(height_matrix):
            line = f"Tree_{i+1:2} |"
            line += "".join(f"{int(val):8d}" for val in row)
            f.write(line + "\n")
        
        # Normalized Matrix
        f.write("\nNormalized Height Matrix:\n\n")
        f.write(header + "\n")
        f.write("-" * (9 + 8 * len(nodes)) + "\n")
        
        for i, row in enumerate(normalized_matrix):
            line = f"Tree_{i+1:2} |"
            line += "".join(f"{val:8.2f}" for val in row)
            f.write(line + "\n")
        
        # Pairwise Differences
        f.write("\nPairwise Differences:\n\n")
        for tree_pair, diff in differences:
            f.write(f"{tree_pair}:\n")
            line = "".join(f"{val:8.2f}" for val in diff)
            f.write(line + "\n\n")

def main():
    base_path = os.path.expanduser("~/1.mahsa.farnia/classificataion_journal")
    input_file = os.path.join(base_path, "weighted_newicks_60.txt")
    output_file = os.path.join(base_path, "height_matrices_60.txt")
    
    with open(input_file, 'r') as f:
        trees = [line.strip() for line in f if line.strip()]
    
    height_matrix, nodes = create_height_matrix(trees)
    normalized_matrix = normalize_matrix(height_matrix)
    differences = calculate_differences(normalized_matrix)
    print_matrices(height_matrix, normalized_matrix, differences, nodes, output_file)

if __name__ == "__main__":
    main()