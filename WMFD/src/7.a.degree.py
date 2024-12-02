import os
import re
import numpy as np

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
    """Get the sequence number for ordering. Returns -1 for 'naive', number for 'seqN'."""
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
                # Read the label after closing parenthesis
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
        
        # Set up parent-child relationships
        node.children = children
        for child in children:
            child.parent = node
            
        return node, pos

    # Clean the Newick string
    newick_str = newick_str.strip()
    if not newick_str.endswith(';'):
        newick_str += ';'
    
    # Parse the tree
    if newick_str.startswith('('):
        root, _ = create_tree_from_string(newick_str, 0)
    else:
        root = Node(newick_str[:-1])
    
    return root

def is_valid_node_name(name):
    """Check if the node name is one we want to include (naive or seqN)."""
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

    # Process direct children
    for child in node.children:
        if not process_child(child):
            # If child is not a valid node, check its children
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
        """Process each node in the tree and calculate its degree."""
        if node.name and is_valid_node_name(node.name):
            base_name = clean_node_name(node.name)
            
            if base_name == 'naive':
                # For naive, count all first-level sequence nodes and unnamed nodes with length
                degrees[base_name] = get_all_first_level_seq_nodes(node)
            else:
                # For sequence nodes, count direct children and unnamed nodes with length
                num_children = count_direct_children(node)
                degrees[base_name] = num_children
        
        # Process all children
        for child in node.children:
            process_node(child)
    
    process_node(root)
    return degrees

def process_newick_file(file_path):
    """Process multiple Newick trees from a file and return their node degrees."""
    all_trees_degrees = []
    all_node_names = set()
    
    with open(file_path, 'r') as f:
        newick_strings = f.readlines()
    
    print(f"Found {len(newick_strings)} lines in the file.")
    
    # First pass: collect all possible node names
    for line in newick_strings:
        if line.strip():
            tree_str = re.sub(r'^.*?:\s*', '', line.strip())
            root = parse_custom_newick(tree_str)
            
            def collect_names(node):
                if node.name and is_valid_node_name(node.name):
                    all_node_names.add(clean_node_name(node.name))
                for child in node.children:
                    collect_names(child)
            
            collect_names(root)
    
    # Sort node names with naive first, then seq1, seq2, etc.
    sorted_node_names = sorted(all_node_names, key=get_node_number)
    print("Node names in order:", sorted_node_names)
    
    # Second pass: calculate degrees
    for i, line in enumerate(newick_strings, 1):
        if line.strip():
            try:
                tree_str = re.sub(r'^.*?:\s*', '', line.strip())
                root = parse_custom_newick(tree_str)
                
                # Calculate degrees
                degrees = calculate_node_degrees(root)
                
                # Ensure all nodes are present in the degrees dictionary
                full_degrees = {node: degrees.get(node, 0) for node in sorted_node_names}
                all_trees_degrees.append(full_degrees)
                
                print(f"Processed tree {i}")
                
            except Exception as e:
                print(f"Error processing tree {i}: {str(e)}")
                continue
    
    if not all_trees_degrees:
        raise Exception("No trees were successfully processed")
    
    return all_trees_degrees, sorted_node_names

def create_degree_matrix(all_trees_degrees, node_names):
    """Create a matrix of node degrees for all trees."""
    num_trees = len(all_trees_degrees)
    num_nodes = len(node_names)
    
    matrix = np.zeros((num_trees, num_nodes), dtype=int)
    
    for i, tree_degrees in enumerate(all_trees_degrees):
        for j, node_name in enumerate(node_names):
            matrix[i, j] = tree_degrees.get(node_name, 0)
    
    return matrix

def min_max_normalize_matrix(matrix):
    """Perform min-max normalization on each column of the matrix."""
    normalized_matrix = np.zeros_like(matrix, dtype=float)
    
    for j in range(matrix.shape[1]):
        column = matrix[:, j]
        col_min = np.min(column)
        col_max = np.max(column)
        
        if col_max - col_min != 0:
            normalized_matrix[:, j] = (column - col_min) / (col_max - col_min)
        else:
            normalized_matrix[:, j] = 0
            
    return normalized_matrix

def calculate_pairwise_differences(normalized_matrix):
    """Calculate differences between pairs of trees using normalized values."""
    num_trees = normalized_matrix.shape[0]
    comparisons = []
    
    for i in range(num_trees):
        for j in range(i + 1, num_trees):  # start from i+1 to avoid duplicates
            diff = np.abs(normalized_matrix[i] - normalized_matrix[j])
            comparisons.append((i+1, j+1, diff))
    
    return comparisons

def save_all_results(matrix, normalized_matrix, comparisons, node_names, output_path):
    """Save all results to a single file."""
    with open(output_path, 'w') as f:
        # Write original matrix
        f.write("ORIGINAL DEGREE MATRIX\n")
        f.write("====================\n")
        header = ['Tree'] + node_names
        f.write('\t'.join(header) + '\n')
        for i in range(matrix.shape[0]):
            row = [f'Tree_{i+1}'] + [str(x) for x in matrix[i, :]]
            f.write('\t'.join(row) + '\n')
        
        # Write normalized matrix
        f.write("\nNORMALIZED MATRIX (MIN-MAX SCALING)\n")
        f.write("=================================\n")
        f.write('\t'.join(header) + '\n')
        for i in range(normalized_matrix.shape[0]):
            row = [f'Tree_{i+1}'] + [f"{x:.4f}" for x in normalized_matrix[i, :]]
            f.write('\t'.join(row) + '\n')
        
        # Write pairwise differences
        f.write("\nPAIRWISE DIFFERENCES (based on normalized values)\n")
        f.write("============================================\n")
        f.write("Tree_i\tTree_j\t" + '\t'.join(node_names) + '\n')
        for i, j, diff in comparisons:
            row = [f"Tree_{i}", f"Tree_{j}"] + [f"{x:.4f}" for x in diff]
            f.write('\t'.join(row) + '\n')

def main():
    # Input and output paths
    input_path = os.path.expanduser('~/1.mahsa.farnia/classificataion_journal/weighted_newicks_60.txt')
    output_path = os.path.expanduser('~/1.mahsa.farnia/classificataion_journal/tree_degrees_results.txt')
    
    try:
        # Process all trees and get degrees
        print("Processing Newick trees...")
        print(f"Reading from: {input_path}")
        
        if not os.path.exists(input_path):
            raise Exception(f"Input file not found: {input_path}")
        
        all_trees_degrees, node_names = process_newick_file(input_path)
        
        # Create degree matrix
        print("\nCreating degree matrix...")
        matrix = create_degree_matrix(all_trees_degrees, node_names)
        
        # Perform min-max normalization
        print("\nPerforming min-max normalization...")
        normalized_matrix = min_max_normalize_matrix(matrix)
        
        # Calculate pairwise differences using normalized values
        print("\nCalculating pairwise differences...")
        comparisons = calculate_pairwise_differences(normalized_matrix)
        
        # Save all results
        print("\nSaving all results...")
        save_all_results(matrix, normalized_matrix, comparisons, node_names, output_path)
        
        print(f"\nAnalysis complete. All results saved to: {output_path}")
        print(f"Successfully processed {len(all_trees_degrees)} trees with {len(node_names)} nodes")
        
    except Exception as e:
        print(f"\nCritical error: {str(e)}")

if __name__ == "__main__":
    main()