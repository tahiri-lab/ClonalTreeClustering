from ete3 import Tree

def extract_bipartitions(newick_str):
    """
    Extracts all non-trivial bipartitions (splits) from a phylogenetic tree 
    given in Newick format, including singleton partitions (leaf taxa).

    Parameters
    ----------
    newick_str : str
        A tree encoded in Newick format, possibly with named internal nodes.

    Returns
    -------
    List[frozenset]
        A list of bipartitions, each represented as a frozenset of taxon labels.
        The list includes singleton partitions and internal bipartitions, and is 
        sorted by partition size and lexicographic order.
    """

    # Parse the Newick string into an ETE tree object
    t = Tree(newick_str, format=1)

    # Set to collect singleton bipartitions (individual leaf taxa)
    singleton_biparts = set()

    def collect_named_taxa(node):
        """
        Recursively collects singleton bipartitions for all named leaf taxa.
        """
        if node.is_leaf() and node.name:
            singleton_biparts.add(frozenset([node.name]))
        for child in node.children:
            collect_named_taxa(child)

    collect_named_taxa(t)

    def label_internal_nodes_as_leaves(node):
        """
        If an internal node has a name, it is treated as a pseudo-leaf
        by creating a child node with the same name, and clearing the parent name.
        This preserves internal node labels as leaves for downstream analyses.
        """
        if node.name and not node.is_leaf():
            node.add_child(name=node.name)
            node.name = ""
        for child in node.children:
            label_internal_nodes_as_leaves(child)

    label_internal_nodes_as_leaves(t)

    # Retrieve the complete set of taxa (leaf names)
    all_taxa = set(t.get_leaf_names())

    # Initialize bipartitions with singleton taxa
    bipartitions = set(singleton_biparts)

    # Traverse the tree in post-order to collect all valid bipartitions
    for node in t.traverse("postorder"):
        if not node.is_leaf():
            # Get the set of descendant taxa under this node
            leaves = set(node.get_leaf_names())
            # Only retain non-trivial bipartitions (excluding full tree and singleton)
            if 1 < len(leaves) < len(all_taxa):
                bipartitions.add(frozenset(leaves))

    # Return bipartitions sorted by size and lexicographic order of taxa
    return [frozenset(b) for b in sorted(bipartitions, key=lambda x: (len(x), sorted(x)))]



# ----------------- Example usage -----------------
if __name__ == "__main__":
    newick = "(((12)2)9,((1,6)7,(3,5,13)10)4,14,(11)15,(16)8)N;"
    biparts = extract_bipartitions(newick)

    for i, b in enumerate(biparts, 1):
        print(f"{i:02d}. {b}")
