from ete3 import Tree

def get_direct_arcs(tree):
    """
    Extracts all direct arcs (parent → child) from a given phylogenetic tree.

    Parameters:
    -----------
    tree : ete3.Tree
        A tree object parsed with ETE3.

    Returns:
    --------
    arcs : set of tuple
        A set of arcs represented as (parent_name, child_name).
    """
    arcs = set()
    for node in tree.traverse():
        for child in node.children:
            arcs.add((node.name, child.name))  # Arc from parent to child
    return arcs


def compute_jaccard_distance_arcs(newick1, newick2):
    """
    Computes the Jaccard distance between two trees based on their direct arcs.

    The Jaccard distance measures dissimilarity between the sets of arcs
    extracted from each tree:
        Jaccard distance = 1 - (|A ∩ B| / |A ∪ B|)

    Parameters:
    -----------
    newick1 : str
        First tree in Newick format.
    
    newick2 : str
        Second tree in Newick format.

    Returns:
    --------
    jaccard_distance : float
        A value in [0, 1] indicating dissimilarity:
        - 0 means identical arc sets
        - 1 means no common arcs
    """
    tree1 = Tree(newick1, format=1)
    tree2 = Tree(newick2, format=1)

    arcs1 = get_direct_arcs(tree1)
    arcs2 = get_direct_arcs(tree2)

    intersection = arcs1.intersection(arcs2)
    union = arcs1.union(arcs2)

    if len(union) == 0:
        return 0.0  # Avoid division by zero in pathological case

    jaccard_distance = 1 - len(intersection) / len(union)

    # Uncomment for debugging:
    # print(f"Common arcs: {intersection}")
    # print(f"Num common arcs: {len(intersection)}")
    # print(f"Num total arcs (union): {len(union)}")

    return jaccard_distance

# Example usage
if __name__ == "__main__":
    newick_str1 = "(((12)2)9,((6,1)7,(3,5,13)10)4,14,(11)15,(16)8)N;"
    newick_str2 = "(((2)9)12,((1,4)7,(10,5,13)6)3,14,(15)11,(16)8)N;"
    newick_str3 = "(((12)9)2,((13,5)1,(3,4,6)7)10,14,(11)15,(16)8)N;"
    newick_str4 = "(((12)2)9,((3,1)7,(13,6,5)4)10,14,(11)15,(8)16)N;"
    newick_str5 = "(((9)2)12,((3,1)7,(13,4,6)5)10,14,(11)15,(8)16)N;"

    distance = compute_jaccard_distance_arcs(newick_str2, newick_str3)
    print(f"\nJaccard distance (on direct arcs): {distance:.4f}")
