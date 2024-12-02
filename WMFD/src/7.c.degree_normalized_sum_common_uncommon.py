import os
import re
import numpy as np
from itertools import combinations

class Node:
    def __init__(self, name=None):
        self.name = name
        self.children = []
        self.parent = None

def clean_node_name(name):
    """Extract the base name before @ symbol."""
    if not name:
        return None
    parts = name.split('@')[0].strip()
    return parts if parts else None

def get_node_number(name):
    """Get the sequence number for ordering."""
    if name == 'naive':
        return -1
    if name.startswith('seq'):
        try:
            return int(name[3:])
        except ValueError:
            return float('inf')
    return float('inf')

def parse_custom_newick(newick_str):
    """Parse the custom Newick format with @ symbols."""
    def create_tree_from_string(s, pos):
        node = Node()
        children = []
        current_token = []
        
        while pos < len(s):
            char = s[pos]
            
            if char == '(':
                child_node, new_pos = create_tree_from_string(s, pos + 1)
                children.append(child_node)
                pos = new_pos
            elif char == ',':
                if current_token:
                    node_name = ''.join(current_token).strip()
                    if node_name:
                        new_node = Node(node_name)
                        children.append(new_node)
                    current_token = []
                pos += 1
            elif char == ')':
                if current_token:
                    node_name = ''.join(current_token).strip()
                    if node_name:
                        new_node = Node(node_name)
                        children.append(new_node)
                    current_token = []
                pos += 1
                while pos < len(s) and s[pos] not in ['(', ')', ',', ';']:
                    current_token.append(s[pos])
                    pos += 1
                if current_token:
                    node.name = ''.join(current_token).strip()
                break
            elif char == ';':
                if current_token:
                    node_name = ''.join(current_token).strip()
                    if node_name:
                        new_node = Node(node_name)
                        children.append(new_node)
                break
            else:
                current_token.append(char)
                pos += 1
        
        node.children = children
        for child in children:
            child.parent = node
            
        return node, pos

    newick_str = newick_str.strip()
    if not newick_str.endswith(';'):
        newick_str += ';'
    
    if newick_str.startswith('('):
        root, _ = create_tree_from_string(newick_str, 0)
    else:
        root = Node(newick_str[:-1])
    
    return root

def is_valid_node_name(name):
    """Check if the node name is valid."""
    if not name:
        return False
    base_name = clean_node_name(name)
    return base_name == 'naive' or base_name.startswith('seq')

def has_branch_length(name):
    """Check if a node has a branch length (contains ':')."""
    return name and ':' in name

def get_all_first_level_seq_nodes(node):
    """Get all sequence nodes and unnamed nodes with branch length at the first level."""
    first_level_nodes = set()
    unnamed_with_length = 0
    
    def process_child(child):
        if child.name:
            if is_valid_node_name(child.name):
                base_name = clean_node_name(child.name)
                if base_name.startswith('seq'):
                    first_level_nodes.add(base_name)
                return True
            elif has_branch_length(child.name):
                nonlocal unnamed_with_length
                unnamed_with_length += 1
                return True
        return False

    for child in node.children:
        if not process_child(child):
            for grandchild in child.children:
                process_child(grandchild)
                
    return len(first_level_nodes) + unnamed_with_length

def count_direct_children(node):
    """Count direct children including valid nodes and unnamed nodes with branch length."""
    count = 0
    for child in node.children:
        if child.name:
            if is_valid_node_name(child.name):
                base_name = clean_node_name(child.name)
                if base_name.startswith('seq'):
                    count += 1
            elif has_branch_length(child.name):
                count += 1
    return count

def calculate_node_degrees(root):
    """Calculate degrees for all named nodes (naive and sequences) in the tree."""
    degrees = {}
    
    def process_node(node):
        if node.name and is_valid_node_name(node.name):
            base_name = clean_node_name(node.name)
            
            if base_name == 'naive':
                degrees[base_name] = get_all_first_level_seq_nodes(node)
            else:
                num_children = count_direct_children(node)
                degrees[base_name] = num_children
        
        for child in node.children:
            process_node(child)
    
    process_node(root)
    return degrees

def min_max_normalize_matrix(matrix):
    """Perform min-max normalization on each column of the matrix."""
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

def analyze_tree_differences():
    # Input and output paths
    input_path = os.path.expanduser('~/1.mahsa.farnia/classificataion_journal/weighted_newicks_60.txt')
    output_path = os.path.expanduser('~/1.mahsa.farnia/classificataion_journal/degree_comparison_normalized_results.txt')
    
    # Read trees
    with open(input_path, 'r') as f:
        trees = [re.sub(r'^.*?:\s*', '', line.strip()) for line in f if line.strip()]
    
    # Extract nodes and calculate degrees
    all_nodes = set()
    tree_nodes = []
    tree_degrees = []
    
    for tree in trees:
        root = parse_custom_newick(tree)
        degrees = calculate_node_degrees(root)
        nodes = set(degrees.keys())
        all_nodes.update(nodes)
        tree_nodes.append(nodes)
        tree_degrees.append(degrees)
    
    # Sort nodes
    sorted_nodes = sorted(all_nodes, key=get_node_number)
    
    # Create degree matrix
    degree_matrix = np.zeros((len(trees), len(sorted_nodes)))
    
    for i, degrees in enumerate(tree_degrees):
        for j, node in enumerate(sorted_nodes):
            degree_matrix[i, j] = degrees.get(node, 0)
    
    # Normalize the matrix
    normalized_matrix = min_max_normalize_matrix(degree_matrix)
    
    # Calculate differences for common and uncommon nodes
    results = []
    for (i, nodes1), (j, nodes2) in combinations(enumerate(tree_nodes), 2):
        # Find common and uncommon nodes
        common_nodes = nodes1 & nodes2
        all_nodes = nodes1 | nodes2
        uncommon_nodes = all_nodes - common_nodes
        
        num_common = len(common_nodes)
        num_uncommon = len(uncommon_nodes)
        
        # Get indices for common and uncommon nodes
        common_indices = [k for k, node in enumerate(sorted_nodes) if node in common_nodes]
        uncommon_indices = [k for k, node in enumerate(sorted_nodes) if node in uncommon_nodes]
        
        # Calculate normalized sums
        common_sum = (np.abs(normalized_matrix[i, common_indices] - 
                           normalized_matrix[j, common_indices]).sum() / 
                     num_common) if num_common > 0 else 0
        
        uncommon_sum = (np.abs(normalized_matrix[i, uncommon_indices] - 
                             normalized_matrix[j, uncommon_indices]).sum() / 
                       num_uncommon) if num_uncommon > 0 else 0
        
        results.append({
            'trees': f"Tree_{i+1}_vs_Tree_{j+1}",
            'common_sum': common_sum,
            'uncommon_sum': uncommon_sum,
            'num_common': num_common,
            'num_uncommon': num_uncommon
        })
    
    # Save results
    with open(output_path, 'w') as f:
        f.write("Normalized Degree Differences (Divided by Node Count):\n")
        f.write("Tree Pair              Common Nodes Sum  Uncommon Nodes Sum  #Common Nodes  #Uncommon Nodes\n")
        f.write("-" * 85 + "\n")
        for result in results:
            f.write(f"{result['trees']:<20} {result['common_sum']:>14.4f} "
                   f"{result['uncommon_sum']:>16.4f} {result['num_common']:>13d} "
                   f"{result['num_uncommon']:>15d}\n")
    
    print(f"Results have been saved to: {output_path}")
    return results

if __name__ == "__main__":
    analyze_tree_differences()
    