import numpy as np
import re
from pathlib import Path
from itertools import combinations
from datetime import datetime

def parse_node_info(s):
    """Extract node name without the @...:... part."""
    if '@' in s:
        return s.split('@')[0]
    return s

def get_all_connections(newick_str):
    """Extract all parent-child connections from Newick string."""
    # Remove the prefix (e.g., "60_1: ") by splitting at the first space
    if ' ' in newick_str:
        newick_str = newick_str.split(' ', 1)[1]
    
    # Remove trailing semicolon and spaces
    newick_str = newick_str.strip().strip(';')
    
    connections = []
    
    def split_top_level(s):
        """Split string by commas at top level only."""
        parts = []
        current = ""
        level = 0
        for c in s:
            if c == '(':
                level += 1
            elif c == ')':
                level -= 1
            elif c == ',' and level == 0 and current:
                parts.append(current)
                current = ""
                continue
            current += c
        if current:  # Add the last part
            parts.append(current)
        return parts

    def process_node(node_str):
        """Process a single node string and return its name."""
        if '@' in node_str:
            return node_str.split('@')[0].strip()
        return node_str.split(':')[0].strip() if ':' in node_str else node_str.strip()

    def extract_node_connections(s, parent=None):
        """Recursively extract all connections from a tree string."""
        local_connections = []
        
        if s.startswith('('):
            # Find matching closing parenthesis
            level = 0
            end_idx = -1
            for i, c in enumerate(s):
                if c == '(':
                    level += 1
                elif c == ')':
                    level -= 1
                    if level == 0:
                        end_idx = i
                        break
            
            inner = s[1:end_idx]
            remainder = s[end_idx + 1:]
            
            # Get current node name
            current_node = process_node(remainder) if remainder else parent
            
            # Process all children
            children = split_top_level(inner)
            for child in children:
                if '(' in child:
                    # This is a subtree
                    local_connections.extend(extract_node_connections(child, current_node))
                else:
                    # This is a leaf node
                    child_name = process_node(child)
                    if current_node and child_name:
                        local_connections.append((current_node, child_name))
            
            if parent and current_node:
                local_connections.append((parent, current_node))
                
        else:
            # This is a leaf node
            node_name = process_node(s)
            if parent and node_name:
                local_connections.append((parent, node_name))
        
        return local_connections

    if newick_str.startswith('('):
        # Process the main content
        connections = extract_node_connections(newick_str, 'naive')
    
    return connections

def create_adjacency_matrix(connections, max_seq_num=28):
    """Create adjacency matrix showing parent-child relationships."""
    nodes = ['naive'] + [f'seq{i}' for i in range(1, max_seq_num + 1)]
    n = len(nodes)
    node_to_idx = {node: i for i, node in enumerate(nodes)}
    matrix = np.zeros((n, n), dtype=int)
    
    for parent, child in connections:
        if parent in node_to_idx and child in node_to_idx:
            parent_idx = node_to_idx[parent]
            child_idx = node_to_idx[child]
            matrix[parent_idx, child_idx] = 1
    
    return matrix

def format_matrix(matrix, string_id):
    """Format matrix as string with headers."""
    nodes = ['naive'] + [f'seq{i}' for i in range(1, 29)]
    output = f"Adjacency Matrix {string_id}:\n"
    output += "    " + " ".join(f"{node:5}" for node in nodes) + "\n"
    
    for i, node in enumerate(nodes):
        output += f"{node:4}" + " ".join(f"{x:5}" for x in matrix[i]) + "\n"
    return output

def count_nodes_in_matrix(matrix):
    """Count the number of nodes (rows that have at least one connection) in a matrix."""
    return np.sum(np.any(matrix != 0, axis=1) | np.any(matrix != 0, axis=0))

def calculate_hamming_distance(matrix1, matrix2):
    """Calculate the Hamming distance between two matrices."""
    return np.sum(matrix1 != matrix2)

def calculate_normalized_distances(matrices):
    """Calculate normalized Hamming distances between all pairs of matrices."""
    num_matrices = len(matrices)
    distances = {}
    
    for (i, j) in combinations(range(num_matrices), 2):
        nodes_i = count_nodes_in_matrix(matrices[i])
        nodes_j = count_nodes_in_matrix(matrices[j])
        normalization_factor = nodes_i + nodes_j - 2
        hamming_dist = calculate_hamming_distance(matrices[i], matrices[j])
        
        normalized_dist = hamming_dist / normalization_factor if normalization_factor > 0 else 0
        
        distances[(i, j)] = {
            'hamming': hamming_dist,
            'normalized': normalized_dist,
            'nodes_i': nodes_i,
            'nodes_j': nodes_j,
            'norm_factor': normalization_factor
        }
    
    return distances

def main():
    # Define file paths
    base_path = Path("/home/local/USHERBROOKE/farm2103/1.mahsa.farnia/classificataion_journal")
    input_file = base_path / "weighted_newicks_60.txt"
    output_file = base_path / "Normalized_HD.txt"

    # Process all the data and write to a single output file
    with open(output_file, 'w') as f:
        # Write header
        f.write("=== Complete Newick String Analysis ===\n")
        f.write(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

        # Read and process Newick strings
        with open(input_file, 'r') as infile:
            newick_strings = [line.strip() for line in infile if line.strip()]

        matrices = []
        string_ids = []

        # Process each Newick string
        for newick in newick_strings:
            string_id = newick.split(':')[0].strip()
            string_ids.append(string_id)
            
            # Process the string
            connections = get_all_connections(newick)
            matrix = create_adjacency_matrix(connections)
            matrices.append(matrix)

        # Calculate distances
        distances = calculate_normalized_distances(matrices)

        # Write all results
        f.write("=== Normalized Hamming Distances ===\n\n")
        
        # Format as a table
        f.write(f"{'Matrix Pair':<20} {'Raw Distance':<15} {'Normalized Distance':<20}\n")
        f.write("-" * 55 + "\n")
        
        for (i, j), info in sorted(distances.items()):
            pair = f"{string_ids[i]}-{string_ids[j]}"
            f.write(f"{pair:<20} {info['hamming']:<15} {info['normalized']:.4f}\n")

        # Write detailed breakdown for verification
        f.write("\n=== Detailed Calculations ===\n\n")
        for (i, j), info in sorted(distances.items()):
            f.write(f"Pair: {string_ids[i]} - {string_ids[j]}\n")
            f.write(f"  Raw Hamming Distance: {info['hamming']}\n")
            f.write(f"  Nodes in {string_ids[i]}: {info['nodes_i']}\n")
            f.write(f"  Nodes in {string_ids[j]}: {info['nodes_j']}\n")
            f.write(f"  Normalization Factor: {info['norm_factor']}\n")
            f.write(f"  Normalized Distance: {info['normalized']:.4f}\n\n")

    print(f"Analysis complete. Results written to: {output_file}")

if __name__ == "__main__":
    main()