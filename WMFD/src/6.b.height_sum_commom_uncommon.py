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

def get_nodes_from_tree(tree):
    heights = get_node_heights(tree)
    return set(heights.keys())

def analyze_tree_differences():
    # Read the input file
    base_path = os.path.expanduser("~/1.mahsa.farnia/classificataion_journal")
    input_file = os.path.join(base_path, "weighted_newicks_60.txt")
    
    with open(input_file, 'r') as f:
        trees = [line.strip() for line in f if line.strip()]
    
    # Get nodes for each tree
    tree_nodes = [get_nodes_from_tree(tree) for tree in trees]
    
    # Calculate height matrices as before
    height_matrix, nodes = create_height_matrix(trees)
    normalized_matrix = normalize_matrix(height_matrix)
    
    # Create result matrix
    n_trees = len(trees)
    results = []
    
    for i in range(n_trees):
        for j in range(i + 1, n_trees):
            # Find common and uncommon nodes
            common_nodes = tree_nodes[i] & tree_nodes[j]
            uncommon_nodes = tree_nodes[i] ^ tree_nodes[j]  # symmetric difference
            
            # Get indices for common and uncommon nodes
            common_indices = [k for k, node in enumerate(nodes) if node in common_nodes]
            uncommon_indices = [k for k, node in enumerate(nodes) if node in uncommon_nodes]
            
            # Calculate differences for this pair of trees
            diff = np.abs(normalized_matrix[i] - normalized_matrix[j])
            
            # Sum the differences for common and uncommon nodes
            common_sum = np.sum(diff[common_indices]) if common_indices else 0
            uncommon_sum = np.sum(diff[uncommon_indices]) if uncommon_indices else 0
            
            results.append({
                'pair': f"Tree_{i+1}-Tree_{j+1}",
                'common_sum': common_sum,
                'uncommon_sum': uncommon_sum
            })
    
    return results

def write_results_to_file(results, output_file):
    with open(output_file, 'w') as f:
        f.write("Results Matrix:\n")
        f.write("-" * 60 + "\n")
        f.write(f"{'Tree Pair':<15} {'Common Nodes Sum':>20} {'Uncommon Nodes Sum':>20}\n")
        f.write("-" * 60 + "\n")
        
        for result in results:
            f.write(f"{result['pair']:<15} {result['common_sum']:20.4f} {result['uncommon_sum']:20.4f}\n")

def main():
    base_path = os.path.expanduser("~/1.mahsa.farnia/classificataion_journal")
    output_file = os.path.join(base_path, "height_sums_common_uncommon.txt")
    
    results = analyze_tree_differences()
    write_results_to_file(results, output_file)
    print(f"Results have been written to: {output_file}")

if __name__ == "__main__":
    main()