# -*- coding: utf-8 -*-
"""
Step 1b) Metric (WMFD) from trees.csv (columns: cluster, tree_id, newick)

Outputs:
- wmfd_matrix.csv: NxN matrix (headers = ids "c.t")
- wmfd_pairs.csv: triples (id_i,id_j,dist) i<j
- labels.csv: ground truth (id, true_cluster)

WMFD = P * WND_uncommon + WND_common + L5 * HD 
P = 1 - |common leaves| / |union leaves| 
WND_* = L1*BL + L2*H + L3*W + L4*D 
HD = 1 - |splits_intersection| / |splits_union|

Notes:
- BL = leaf branch length
- H = root→leaf distance
- W = 1.0 (default weight)
- D = 0.0 for leaves (keeps location for future compatibility)
- We cache features and splits for performance.
"""

import csv
import argparse
from typing import Dict, Set, List, Tuple
from ete3 import Tree

def read_trees_csv(in_csv: str):
    ids: List[str] = []
    clusters: List[int] = []
    newicks: List[str] = []
    with open(in_csv, newline="", encoding="utf-8") as f:
        r = csv.DictReader(f)
        for row in r:
            c = int(row["cluster"])
            t = int(row["tree_id"])
            nwk = row["newick"].strip()
            ids.append(f"{c}.{t}")
            clusters.append(c)
            newicks.append(nwk)
    return ids, clusters, newicks

def precompute_features(nwk: str):
    """Retourne (tree, leaf_BL, leaf_H, leaf_W, leaf_D, leaf_set, splits_set)."""
    t = Tree(nwk, format=1)
    root = t.get_tree_root()

    # features leaves
    leaf_BL: Dict[str, float] = {}
    leaf_H: Dict[str, float] = {}
    leaf_W: Dict[str, float] = {}
    leaf_D: Dict[str, float] = {}
    for leaf in t.iter_leaves():
        name = str(leaf.name)
        leaf_BL[name] = float(leaf.dist or 0.0)
        leaf_H[name]  = float(root.get_distance(leaf) or 0.0)
        leaf_W[name]  = 1.0
        leaf_D[name]  = 0.0  # placeholder (0 for leaves)

    leaf_set: Set[str] = set(leaf_BL.keys())

    # splits (set of internal clades as sets of sheets)
    splits: Set[frozenset] = set()
    for node in t.traverse("postorder"):
        if node.is_leaf():
            continue
        clade = frozenset(l.name for l in node.iter_leaves())
        if 1 < len(clade) < len(leaf_set):
            splits.add(clade)

    return t, leaf_BL, leaf_H, leaf_W, leaf_D, leaf_set, splits

def wmfd_pair(
    A, B,
    L1: float, L2: float, L3: float, L4: float, L5: float
) -> float:
    """
    A, B: tuples retournés par precompute_features(...)
    """
    (_t1, bl1, h1, w1, d1, Ls1, S1) = A
    (_t2, bl2, h2, w2, d2, Ls2, S2) = B

    unionL = Ls1 | Ls2
    interL = Ls1 & Ls2
    TN = len(unionL)
    CN = len(interL)
    P = 1.0 - (CN / TN) if TN > 0 else 0.0

    # agrégats common / uncommon
    com_BL = com_H = com_W = com_D = 0.0
    unc_BL = unc_H = unc_W = unc_D = 0.0

    for leaf in unionL:
        # missing values => 0.0
        diff_BL = abs(bl1.get(leaf, 0.0) - bl2.get(leaf, 0.0))
        diff_H  = abs(h1.get(leaf, 0.0)  - h2.get(leaf, 0.0))
        diff_W  = abs(w1.get(leaf, 0.0)  - w2.get(leaf, 0.0))
        diff_D  = abs(d1.get(leaf, 0.0)  - d2.get(leaf, 0.0))
        if leaf in interL:
            com_BL += diff_BL; com_H += diff_H; com_W += diff_W; com_D += diff_D
        else:
            unc_BL += diff_BL; unc_H += diff_H; unc_W += diff_W; unc_D += diff_D

    BL_common = (com_BL / CN) if CN > 0 else 0.0
    H_common  = (com_H  / CN) if CN > 0 else 0.0
    W_common  = (com_W  / CN) if CN > 0 else 0.0
    D_common  = (com_D  / CN) if CN > 0 else 0.0

    nb_unc = TN - CN
    BL_uncommon = (unc_BL / nb_unc) if nb_unc > 0 else 0.0
    H_uncommon  = (unc_H  / nb_unc) if nb_unc > 0 else 0.0
    W_uncommon  = (unc_W  / nb_unc) if nb_unc > 0 else 0.0
    D_uncommon  = (unc_D  / nb_unc) if nb_unc > 0 else 0.0

    WND_common   = L1*BL_common   + L2*H_common   + L3*W_common   + L4*D_common
    WND_uncommon = L1*BL_uncommon + L2*H_uncommon + L3*W_uncommon + L4*D_uncommon

    # HD on splits
    unionS = S1 | S2
    interS = S1 & S2
    HD = 1.0 - (len(interS)/len(unionS)) if len(unionS) > 0 else 0.0

    return float(P * WND_uncommon + WND_common + L5 * HD)

def write_outputs(ids, clusters, D, out_matrix, out_pairs, out_labels):
    n = len(ids)
    # matrix
    with open(out_matrix, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["id"] + ids)
        for i in range(n):
            w.writerow([ids[i]] + [("{:.8f}".format(D[i][j])) for j in range(n)])
    # pairs
    with open(out_pairs, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["id_i","id_j","dist"])
        for i in range(n):
            for j in range(i+1, n):
                w.writerow([ids[i], ids[j], "{:.8f}".format(D[i][j])])
    # labels (field truth)
    with open(out_labels, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["id","true_cluster"])
        for i, lab in enumerate(clusters):
            w.writerow([ids[i], lab])

def main():
    ap = argparse.ArgumentParser(description="WMFD pairwise distances from trees.csv")
    ap.add_argument("--in_csv", type=str, required=True, help="trees.csv (cluster,tree_id,newick)")
    ap.add_argument("--out_matrix", type=str, default="wmfd_matrix.csv")
    ap.add_argument("--out_pairs",  type=str, default="wmfd_pairs.csv")
    ap.add_argument("--out_labels", type=str, default="labels.csv")
    # weightings
    ap.add_argument("--L1", type=float, default=0.30)
    ap.add_argument("--L2", type=float, default=0.20)
    ap.add_argument("--L3", type=float, default=0.25)
    ap.add_argument("--L4", type=float, default=0.15)
    ap.add_argument("--L5", type=float, default=0.10)
    args = ap.parse_args()

    ids, clusters, newicks = read_trees_csv(args.in_csv)

    # pre-calculation of the features
    feats = [precompute_features(nwk) for nwk in newicks]

    n = len(ids)
    D = [[0.0]*n for _ in range(n)]
    for i in range(n):
        if (i % 5) == 0:
            print(f"[WMFD] progress {i}/{n}")
        for j in range(i+1, n):
            d = wmfd_pair(feats[i], feats[j], args.L1, args.L2, args.L3, args.L4, args.L5)
            D[i][j] = D[j][i] = d
    print(f"[WMFD] progress {n}/{n} ✓")

    write_outputs(ids, clusters, D, args.out_matrix, args.out_pairs, args.out_labels)
    print(f"OK. Wrote:\n- {args.out_matrix}\n- {args.out_pairs}\n- {args.out_labels}")

if __name__ == "__main__":
    main()
