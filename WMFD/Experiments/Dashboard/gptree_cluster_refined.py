# -*- coding: utf-8 -*-
"""
gptree_cluster_refined.py — (1a) generator aligned with loops (k, L, Ngen, plevel)

- K=1 fast-path (no overlap constraint, never failure)
- NO blocking timeout: no TimeoutError; we have a maximum number of tries per tree,
  then fallback to the best candidate (closest to Plevel)
- Progressive internal tolerance 
- Avoids reusing the same species tree topology between clusters (k>=2)

CSV output: cluster,tree_id,newick
"""


import argparse
import csv
import random
import sys
from typing import List, Tuple, Set, FrozenSet

from ete3 import Tree
import asymmetree.treeevolve as te
from asymmetree.tools.PhyloTreeTools import to_newick

# ---- internal tolerance (not exposed) ----
EPS_INTERNAL = 0.04  

# ---- validations -------------------------------------------------
def validate_args(k: int, L: int, Ngen: int, plevel: float):
    if not (1 <= k <= 100):
        raise ValueError("k must be in [1,100].")
    if not (5 <= L <= 500):
        raise ValueError("L must be in [5,500].")
    if not (3 <= Ngen <= 1000):
        raise ValueError("Ngen must be in [3,1000].")
    
    if not (0.3 <= plevel <= 0.7):
        raise ValueError("plevel (overlap) must be in [0.3,0.7].")

# ---- tools ------------------------------------------------------
def rename_leaves(ete_tree: Tree) -> Tree:
    for leaf in ete_tree.iter_leaves():
        if not leaf.name.startswith("L"):
            leaf.name = f"L{leaf.name}"
    return ete_tree

def gptree_speciestree(L: int):
    return te.species_tree_n(L)

def gptree_genetree(S, hgt_rate=0.2, loss_rate=0.2, replace_prob=0.9) -> Tree:
    sim = te.GeneTreeSimulator(S)
    T = sim.simulate(hgt_rate=hgt_rate, loss_rate=loss_rate, replace_prob=replace_prob)
    og = te.prune_losses(T)
    return rename_leaves(Tree(to_newick(og, reconc=False), format=1))

def jaccard_leaf_overlap(t1: Tree, t2: Tree) -> float:
    s1 = set(t1.get_leaf_names()); s2 = set(t2.get_leaf_names())
    if not s1 and not s2:
        return 0.0
    inter = len(s1 & s2)
    uni   = len(s1 | s2)
    return inter / float(uni) if uni > 0 else 0.0

def avg_overlap_with_cluster(candidate: Tree, cluster: List[Tree]) -> float:
    if not cluster:
        return 1.0
    vals = [jaccard_leaf_overlap(candidate, t) for t in cluster]
    return sum(vals) / len(vals)

# ---- topological signatures (to avoid duplicates of S) ------
def tree_topology_signature(ete_tree: Tree) -> Tuple[FrozenSet[str], FrozenSet[FrozenSet[str]]]:
    leafset = frozenset(ete_tree.get_leaf_names())
    splits: Set[FrozenSet[str]] = set()
    for node in ete_tree.traverse("postorder"):
        if node.is_leaf():
            continue
        clade = frozenset(leaf.name for leaf in node.iter_leaves())
        if 1 < len(clade) < len(leafset):
            splits.add(clade)
    return (leafset, frozenset(splits))

def species_topology_signature(S) -> Tuple[FrozenSet[str], FrozenSet[FrozenSet[str]]]:
    nwk = to_newick(S, reconc=False)
    ete_t = Tree(nwk, format=1)
    return tree_topology_signature(ete_t)

# ---- building a cluster without blocking timeout ------------
def build_cluster(
    species_tree,
    Ngen: int,
    plevel: float,
    max_tries_per_tree: int,
    hgt_rate: float,
    loss_rate: float,
    replace_prob: float,
) -> List[Tree]:
    """
    Construit un cluster de Ngen arbres.
    
    """
    cluster: List[Tree] = []

    # 1st arbitrary tree
    cluster.append(gptree_genetree(species_tree, hgt_rate, loss_rate, replace_prob))

    while len(cluster) < Ngen:
        accepted = False
        best_tree = None
        best_diff = 1e9

        for attempt in range(1, max_tries_per_tree + 1):
            cand = gptree_genetree(species_tree, hgt_rate, loss_rate, replace_prob)
            ov = avg_overlap_with_cluster(cand, cluster)

            # gradually widens the acceptance window
            widen = 0.01 * (attempt // 100)  # +0.01 every 100 tries
            eps_eff = EPS_INTERNAL + widen

            if (plevel - eps_eff) <= ov <= (plevel + eps_eff):
                cluster.append(cand)
                accepted = True
                break

            diff = abs(ov - plevel)
            if diff < best_diff:
                best_diff = diff
                best_tree = cand

        if not accepted:
            # fallback : we are still moving forward with the best candidate we met
            cluster.append(best_tree if best_tree is not None else gptree_genetree(species_tree, hgt_rate, loss_rate, replace_prob))

    return cluster

# ---- main -------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description="Generate clusters of phylogenetic trees with target overlap (1a).")
    p.add_argument("--k", type=int, required=True, help="number of clusters")
    p.add_argument("--L", type=int, required=True, help="number of leaves for species tree")
    p.add_argument("--Ngen", type=int, required=True, help="trees per cluster")
    p.add_argument("--plevel", type=float, required=True, help="target average leaf-overlap (0.3–0.7)")

    p.add_argument("--seed", type=int, default=0, help="random seed")
    p.add_argument("--out", type=str, required=True, help="output CSV (cluster,tree_id,newick)")

    # kept for compatibility, but ignored on the blocking logic side
    p.add_argument("--timeout_s", type=int, default=0, help="(ignored) kept for backward-compat; no hard timeout")
    p.add_argument("--max_tries_per_tree", type=int, default=2000, help="max attempts before fallback")

    p.add_argument("--hgt", type=float, default=0.2, help="HGT rate for simulator")
    p.add_argument("--loss", type=float, default=0.2, help="loss rate for simulator")
    p.add_argument("--replace_prob", type=float, default=0.9, help="replacement prob for simulator")

    args = p.parse_args()

    try:
        validate_args(args.k, args.L, args.Ngen, args.plevel)
    except ValueError as e:
        print(e)
        sys.exit(1)

    # seed Python
    random.seed(args.seed or 0)

    # K=1 : fast-path simple, never failure, no overlap constraint
    if args.k == 1:
        S = gptree_speciestree(args.L)
        trees = [gptree_genetree(S, args.hgt, args.loss, args.replace_prob) for _ in range(args.Ngen)]
        with open(args.out, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["cluster", "tree_id", "newick"])
            for i, t in enumerate(trees, 1):
                w.writerow([1, i, t.write(format=1)])
        print(f"[K=1 fast-path] Wrote {args.out} (N={args.Ngen})")
        return

    # K >= 2: ensures different species topologies between clusters
    seen_specs: Set[Tuple[FrozenSet[str], FrozenSet[FrozenSet[str]]]] = set()

    with open(args.out, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["cluster", "tree_id", "newick"])

        for c in range(1, args.k + 1):
            # draws a new species topology not yet seen
            while True:
                S = gptree_speciestree(args.L)
                sig = species_topology_signature(S)
                if sig not in seen_specs:
                    seen_specs.add(sig)
                    break

            # little 'auto-relax' for tall trees (more similar trees => easier)
            hgt = args.hgt
            loss = args.loss
            rprob = args.replace_prob
            if args.plevel >= 0.55:
                if hgt > 0.15: hgt = 0.15
                if loss > 0.15: loss = 0.15
                if rprob < 0.95: rprob = 0.95

            trees_c = build_cluster(
                species_tree=S,
                Ngen=args.Ngen,
                plevel=args.plevel,
                max_tries_per_tree=args.max_tries_per_tree,
                hgt_rate=hgt,
                loss_rate=loss,
                replace_prob=rprob,
            )

            for i, t in enumerate(trees_c, 1):
                w.writerow([c, i, t.write(format=1)])

    print(f"[OK] Wrote {args.out} with {args.k}×{args.Ngen} trees (L={args.L}, p={args.plevel})")

if __name__ == "__main__":
    main()
