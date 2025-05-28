from ete3 import Tree
from collections import defaultdict, Counter
from average_Jaccard_distance import compute_average_jaccard_distance
from frequent_arc import extract_frequent_arcs
from Reference_tab import get_direct_arcs

def group_arcs_by_child(newick_list, frequent_arcs):
    """
    Groups non-frequent arcs by their child node and counts how many trees each arc appears in.
    Only includes arcs whose child is NOT already assigned through a frequent arc.

    Parameters:
    -----------
    newick_list : list of str
        A list of phylogenetic trees in Newick format.

    frequent_arcs : set of tuple
        A set of arcs (parent_name, child_name) considered frequent and thus excluded 
        from this analysis.

    Returns:
    --------
    grouped_arcs : dict
        A dictionary mapping each unassigned child node to a Counter that tracks 
        how often each (parent, child) arc appears across the tree set:
        {
            child_name: {
                (parent, child): count_in_trees,
                ...
            },
            ...
        }
    """
    arc_occurrences = defaultdict(Counter)

    # Identify all children already covered by frequent arcs
    assigned_children = {child for _, child in frequent_arcs}

    for newick in newick_list:
        tree = Tree(newick, format=1)  # Read each tree
        arcs = get_direct_arcs(tree)  # Get all direct arcs from the tree
        non_frequent_arcs = arcs - frequent_arcs  # Filter out frequent arcs

        # Count non-frequent arcs by their child
        for arc in non_frequent_arcs:
            parent, child = arc
            if child not in assigned_children:
                arc_occurrences[child][arc] += 1

    return dict(arc_occurrences)


# Example usage
if __name__ == "__main__":
    trees = [
        "(((12)2)9,((1,4)7,(5,6,13)3)10,14,(11)15,(16)8)N;",
        "(((2)9)12,((3,1)7,(13,5,6)4)10,14,(15)11,(8)16)N;",
        "(((12)2)9,((13,4)7,(5,6,1)3)10,14,(15)11,(16)8)N;",
        "(((9)2)12,((1,3)7,(4,5,13)6)10,14,(11)15,(8)16)N;",
        "(((12)9)2,((6,1)7,(13,4,5)3)10,14,(11)15,(16)8)N;"
    ]

    # Step 1: Compute the frequency threshold using average Jaccard distance (L * n rule)
    threshold = compute_average_jaccard_distance(trees, verbose=False)

    # Step 2: Extract arcs that appear frequently enough across trees
    frequent_arcs = extract_frequent_arcs(trees, threshold, verbose=False)

    # Step 3: Group non-frequent arcs by child, filtering out those assigned in frequent arcs
    grouped_arcs = group_arcs_by_child(trees, frequent_arcs)

    # Step 4: Display grouped arcs with counts
    for child, arc_counts in grouped_arcs.items():
        print(f"\nChild: {child}")
        for arc, count in arc_counts.items():
            print(f"  Arc {arc} appears in {count} tree(s)")
