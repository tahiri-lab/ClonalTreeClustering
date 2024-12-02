import numpy as np
import re
from pathlib import Path

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

# Process input file
input_file = Path("/home/local/USHERBROOKE/farm2103/1.mahsa.farnia/classificataion_journal/weighted_newicks_60.txt")
output_file = Path("/home/local/USHERBROOKE/farm2103/1.mahsa.farnia/classificataion_journal/adjacency_matrices.txt")

# Process all strings
with open(input_file, 'r') as f:
    newick_strings = [line.strip() for line in f if line.strip()]

matrices = []
text_output = ""

for i, newick in enumerate(newick_strings, 1):
    # Extract the identifier (60_1, 60_2, etc.) from the start of the string
    string_id = newick.split(':')[0].strip()
    
    connections = get_all_connections(newick)
    matrix = create_adjacency_matrix(connections)
    matrices.append(matrix)
    
    if i > 1:
        text_output += "\n\n"
    text_output += format_matrix(matrix, string_id)

# Save outputs
with open(output_file, 'w') as f:
    f.write(text_output)

matrices_array = np.array(matrices)
np.save(output_file.parent / "adjacency_matrices.npy", matrices_array)