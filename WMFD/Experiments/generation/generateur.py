# -*- coding: utf-8 -*-


import argparse
import sys
import subprocess
import random
from ete3 import Tree
import asymmetree.treeevolve as te
from asymmetree.tools.PhyloTreeTools import to_newick

def install_packages():
    package = ['pandas', 'ete3', 'PyQt5', 'asymmetree']
    subprocess.check_call([sys.executable, '-m', 'pip', 'install'] + package, stdout=subprocess.DEVNULL)

def validate_args(k, L, Ngen, plevel):
    if not (1 <= k <= 100) or not (5 <= L <= 500) or not (3 <= Ngen <= 500):
        raise ValueError("Invalid parameter values. Please check the range and constraints.")
    if not (0.2 <= plevel <= 0.7):
        raise ValueError("Overlap level (plevel) must be between 0.2 and 0.7.")

# Function to rename leaves of a tree with "L" prefix
def rename_leaves(tree):
    for leaf in tree.iter_leaves():  # iterate over leaves only
        leaf.name = f"L{leaf.name}"
    return tree

def gptree_speciestree(L):
    return te.species_tree_n(L)

# Function to generate one gene tree and rename leaves
def gptree_genetree(S1, hgt_rate=0.2, loss_rate=0.2, replace_prob=0.9):
    tree_simulator = te.GeneTreeSimulator(S1)
    T1 = tree_simulator.simulate(hgt_rate=hgt_rate, loss_rate=loss_rate, replace_prob=replace_prob)
    ogt = te.prune_losses(T1)
    Gn_Tree = Tree(to_newick(ogt, reconc=False), format=1)
    return rename_leaves(Gn_Tree)

def calculate_overlap(tree1, tree2):
    # Jaccard over leaf sets
    leaves1 = set(tree1.get_leaf_names())
    leaves2 = set(tree2.get_leaf_names())
    common = leaves1 & leaves2
    union = leaves1 | leaves2
    return len(common) / len(union) if union else 0.0

def gptree_cluster_gene(sptree, Ngen, plevel, tol=0.03, max_tries=20000):
    cluster_dataset = [gptree_genetree(sptree)]
    #print("Now we have 1 tree")
    tries = 0
    while len(cluster_dataset) < Ngen and tries < max_tries:
        tries += 1
        gene_tree_next = gptree_genetree(sptree)
        overlap_levels = [calculate_overlap(gene_tree_next, gt) for gt in cluster_dataset]
        average_overlap = sum(overlap_levels) / len(overlap_levels)
        if plevel - tol <= average_overlap <= plevel + tol:
            cluster_dataset.append(gene_tree_next)
            #print(f"Now we have {len(cluster_dataset)} trees")
    if len(cluster_dataset) < Ngen:
        raise RuntimeError("N'a pas pu atteindre Ngen avec les paramètres/tolérance donnés.")
    return cluster_dataset

# Topology signature helpers
def tree_topology_signature(ete_tree: Tree):
    """
    Return a hashable topology signature:
    - leafset: frozenset of leaf names
    - splits: frozenset of frozensets, each being the set of leaf names under an internal node
    Branch lengths are ignored.
    """
    leafset = frozenset(ete_tree.get_leaf_names())
    splits = set()
    for node in ete_tree.traverse("postorder"):
        if not node.is_leaf():
            clade = frozenset(leaf.name for leaf in node.iter_leaves())
            if 1 < len(clade) < len(leafset):
                splits.add(clade)
    return (leafset, frozenset(splits))

def species_topology_signature(species_tree_obj):
    """
    Convert the asymmetree species tree to Newick, parse with ETE, and build the signature.
    """
    nwk = to_newick(species_tree_obj, reconc=False)
    ete_t = Tree(nwk, format=1)
    return tree_topology_signature(ete_t)

def main_generateur():
    parser = argparse.ArgumentParser(description="Generate clusters of phylogenetic trees with specified overlap.")
    parser.add_argument("k", type=int, help="Number of clusters (1-100)")
    parser.add_argument("L", type=int, help="Number of leaves (5-499)")
    parser.add_argument("Ngen", type=int, help="Number of trees per cluster (3-500)")
    parser.add_argument("plevel", type=float, help="Average overlap between trees (0.2-0.7)")
    args = parser.parse_args()

    try:
        validate_args(args.k, args.L, args.Ngen, args.plevel)
    except ValueError as e:
        print(e)
        sys.exit(1)

    # Keep track of seen species-tree topologies
    seen_species_topologies = set()
    tree_data: dict[tuple[int, int], Tree] = {}

    for k in range(1, args.k + 1):
        #print(f"For cluster {k}")

        # Cherche une espèce avec topologie nouvelle (avec garde-fou)
        MAX_TRIES_SPECIES = 2000
        for _ in range(MAX_TRIES_SPECIES):
            species_tree_k = gptree_speciestree(args.L)  # <-- FIX: args.L
            sig = species_topology_signature(species_tree_k)
            if sig not in seen_species_topologies:
                seen_species_topologies.add(sig)
                break
        else:
            raise RuntimeError("Impossible de trouver une nouvelle topologie d'espèce.")

        # Génère les arbres de gènes pour ce cluster
        cluster_k = gptree_cluster_gene(species_tree_k, args.Ngen, args.plevel)

        # Stocke chaque arbre dans la structure
        for N, tree in enumerate(cluster_k, start=1):
            #tree_data[(k, N)] = tree
            tree_data[(k, N)] = tree.write()

        

    print(f"All done. In-memory structure 'tree_data' populated with {len(tree_data)} trees.")
    print(tree_data)

if __name__ == "__main__":
    install_packages()
    main_generateur()
