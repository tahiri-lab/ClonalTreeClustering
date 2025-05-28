from ete3 import Tree
from collections import Counter
from average_Jaccard_distance import *
from jaccard_distance import *

def extract_frequent_arcs(newick_list, frequency_threshold, verbose=True):
    """
    Extract arcs (edges from parent to child nodes) that appear in at least
    `frequency_threshold` number of trees.

    Parameters:
    newick_list (list of str): List of trees in Newick format.
    frequency_threshold (int): Minimum number of trees an arc must appear in to be considered frequent.
    verbose (bool): If True, prints diagnostic information about extracted arcs.

    Returns:
    set of tuple: A set of arcs (parent_name, child_name) that meet the frequency threshold.
    """
    n = len(newick_list)
    arc_counter = Counter()

    for newick in newick_list:
        tree = Tree(newick, format=8)
        # Extract all parent → child arcs in the tree
        arcs = {(node.name, child.name) for node in tree.traverse() for child in node.children}
        arc_counter.update(arcs)

    # Select arcs that appear in at least `frequency_threshold` trees
    frequent_arcs = {arc for arc, count in arc_counter.items() if count >= frequency_threshold}

    if verbose:
        print(f"\n[Result] Total unique arcs: {len(arc_counter)}")
        print(f"[Result] Arcs appearing in ≥ {frequency_threshold} trees: {len(frequent_arcs)}\n")
        for arc in sorted(frequent_arcs):
            print(f"  {arc[0]} → {arc[1]}  (count: {arc_counter[arc]})")

    return frequent_arcs


# Example usage
if __name__ == "__main__":
    # A list of phylogenetic trees in Newick format
    trees = [
        "(((12)2)9,((6,1)7,(3,5,13)10)4,14,(11)15,(16)8)N;",
        "(((2)9)12,((1,4)7,(10,5,13)6)3,14,(15)11,(16)8)N;",
        "(((12)9)2,((13,5)1,(3,4,6)7)10,14,(11)15,(16)8)N;",
        "(((12)2)9,((3,1)7,(13,6,5)4)10,14,(11)15,(8)16)N;",
        "(((9)2)12,((3,1)7,(13,4,6)5)10,14,(11)15,(8)16)N;"
    ]

    # Step 1: Compute the average Jaccard distance threshold (L * n rule)
    threshold = compute_average_jaccard_distance(trees)

    # Step 2: Extract arcs appearing in at least `threshold` trees
    frequent_arcs = extract_frequent_arcs(trees, threshold, verbose=True)

    # Display the set of frequent arcs
    print(frequent_arcs)
