import os

def read_newick_trees(file_path):
    """
    Read multiple Newick trees from a text file.
    Returns a list of Newick strings.
    """
    trees = []
    try:
        with open(file_path, 'r') as file:
            for line in file:
                line = line.strip()
                if line:  # Skip empty lines
                    trees.append(line)
    except FileNotFoundError:
        raise FileNotFoundError(f"File not found: {file_path}")
    
    return trees

def extract_sequence_names(newick_str):
    """
    Extract sequence names from a Newick string.
    Returns a set of sequence names (without weights and other annotations).
    Including 'naive' node.
    """
    sequence_names = set()
    parts = newick_str.strip(';').replace('(', ',').replace(')', ',').split(',')
    
    for part in parts:
        part = part.strip()
        if part:
            if '@' in part:
                seq_name = part.split('@')[0]
                sequence_names.add(seq_name)  # Include all nodes, including 'naive'
            elif part == 'naive':
                sequence_names.add(part)
    
    return sequence_names

def calculate_penalty(tree1_str, tree2_str):
    """
    Calculate penalty between two trees based on common and total unique nodes.
    """
    seq1 = extract_sequence_names(tree1_str)
    seq2 = extract_sequence_names(tree2_str)
    
    common_nodes = len(seq1.intersection(seq2))
    all_unique_nodes = len(seq1.union(seq2))
    
    penalty = 1 - (common_nodes / all_unique_nodes)
    return penalty

def main():
    base_path = os.path.expanduser("~/1.mahsa.farnia/classificataion_journal")
    input_file = os.path.join(base_path, "weighted_newicks_60.txt")
    output_file = os.path.join(base_path, "penalties.txt")
    
    newick_strings = read_newick_trees(input_file)
    n_trees = len(newick_strings)
    
    # Calculate penalties and write to file
    with open(output_file, 'w') as f:
        for i in range(n_trees):
            for j in range(i+1, n_trees):
                penalty = calculate_penalty(newick_strings[i], newick_strings[j])
                f.write(f"Penalty(Tree_{i+1}, Tree_{j+1})= {penalty:.4f}\n")
    
    print(f"Results have been saved to: {output_file}")

if __name__ == "__main__":
    main()