# -*- coding: utf-8 -*-
# wmfd_all_from_tree_data.py
# Input : either
#   - dict {(k,N): "<newick string>"}
#   - list ["<newick string>", ...]
# Output: {"D": np.ndarray, "labels": [...]}   (labels style configurable)

from __future__ import annotations
from typing import Dict, Tuple, List, Set, FrozenSet, Union
import numpy as np
from ete3 import Tree

# ---------- extraction depuis n'importe quel newick ----------

def _leaf_features_from_newick(nwk: str):
    """
    Parse un Newick et extrait:
      - W: poids des feuilles (déduit de 'name@weight' si présent, sinon 1.0)
      - BL: longueurs de branches terminales
      - H: hauteurs (distance racine -> feuille)
      - D: degré du parent de la feuille (nb d'enfants du noeud parent)
    """
    t = Tree(nwk, format=1)
    root = t.get_tree_root()

    W, BL, H, D = {}, {}, {}, {}
    for lf in t.iter_leaves():
        nm = str(lf.name) if lf.name is not None else ""
        wt = 1.0
        if "@" in nm:       # format "seq@0.7"
            base, wstr = nm.split("@", 1)
            nm = base.strip()
            try:
                wt = float(wstr)
            except ValueError:
                wt = 1.0

        BL[nm] = float(lf.dist or 0.0)
        H[nm]  = float(root.get_distance(lf) or 0.0)
        W[nm]  = float(wt)
        D[nm]  = float(len(lf.up.children) if lf.up else 1.0)

    return W, BL, H, D

def _collect_all_seqs(nwicks: List[str]) -> List[str]:
    S: Set[str] = set()
    for s in nwicks:
        W, BL, H, D = _leaf_features_from_newick(s)
        S.update(W.keys())
    # tri stable: 'naive' d'abord, puis L<num> / seq<num> si possible
    def key(seq):
        if seq == "naive": return (-2, 0, "")
        if seq.startswith("L"):
            try: return (-1, int(seq[1:]), "")
            except: return (-1, 10**9, seq)
        if seq.startswith("seq"):
            try: return (0, int(seq[3:]), "")
            except: return (0, 10**9, seq)
        return (1, 0, seq)
    return sorted(S, key=key)

def _build_and_normalize(nwicks: List[str], seqs: List[str]):
    """
    Construit 4 matrices (W, BL, H, D) taille (nTrees x nSeqs),
    puis normalise PAR COLONNE (min-max) chaque matrice.
    """
    nT, nS = len(nwicks), len(seqs)
    W  = np.zeros((nT, nS), dtype=float)
    BL = np.zeros((nT, nS), dtype=float)
    H  = np.zeros((nT, nS), dtype=float)
    DG = np.zeros((nT, nS), dtype=float)

    for i, s in enumerate(nwicks):
        wi, bli, hi, di = _leaf_features_from_newick(s)
        for j, name in enumerate(seqs):
            W[i, j]  = wi.get(name, 0.0)
            BL[i, j] = bli.get(name, 0.0)
            H[i, j]  = hi.get(name, 0.0)
            DG[i, j] = di.get(name, 0.0)

    def minmax(X):
        Y = np.zeros_like(X, dtype=float)
        for j in range(X.shape[1]):
            col = X[:, j]
            mn, mx = col.min(), col.max()
            Y[:, j] = 0.0 if mx == mn else (col - mn) / (mx - mn)
        return Y

    return minmax(W), minmax(BL), minmax(H), minmax(DG)

# ---------- commun/uncommon, penalty, topologie ----------

def _topology_splits(nwk: str) -> Set[FrozenSet[str]]:
    t = Tree(nwk, format=1)
    leaves = set(t.get_leaf_names())
    splits: Set[FrozenSet[str]] = set()
    for nd in t.traverse("postorder"):
        if nd.is_leaf(): continue
        clade = frozenset(x.name for x in nd.iter_leaves())
        if 1 < len(clade) < len(leaves):
            splits.add(clade)
    return splits

# ---------- pré-calcul pour accélérer les paires ----------

def _precompute_per_tree(nwicks: List[str]):
    """
    Pour chaque arbre:
      - nodes: set des noms de feuilles
      - splits: ensemble de bipartitions
    (plus tard, les matrices normalisées W/BL/H/D sont déjà calculées par ailleurs)
    """
    nodes_list: List[Set[str]] = []
    splits_list: List[Set[FrozenSet[str]]] = []
    for s in nwicks:
        W, BL, H, D = _leaf_features_from_newick(s)
        nodes_list.append(set(W.keys()))
        splits_list.append(_topology_splits(s))
    return nodes_list, splits_list

def _common_uncommon_indices(seqs: List[str], s1: Set[str], s2: Set[str]):
    C = s1 & s2
    U = s1 | s2
    Ic = [j for j, name in enumerate(seqs) if name in C]
    Iu = [j for j, name in enumerate(seqs) if (name in U and name not in C)]
    return Ic, Iu, C, U

def _sum_abs_diff_rows(X: np.ndarray, i: int, j: int, idxs: List[int], divide_by_count=True):
    if not idxs: return 0.0
    diff = np.abs(X[i, idxs] - X[j, idxs]).sum()
    return float(diff / len(idxs)) if divide_by_count else float(diff)

def _penalty(C: Set[str], U: Set[str]) -> float:
    return 0.0 if len(U) == 0 else (1.0 - (len(C) / len(U)))

# ---------- API principale ----------

def compute_wmfd_all(
    tree_data: Union[Dict[Tuple[int, int], str], List[str]],
    lambdas=(0.15, 0.15, 0.25, 0.35, 0.10),  # (BL, W, Degree, Height, Topo)
    progress: bool = True,
    label_style: str = "auto",  # "auto" | "tuple" | "string" | "index" | "none"
):
    """
    Calcule UNE matrice WMFD pour tous les arbres (avec ou sans '@poids').
    'tree_data' peut être:
      - dict {(k,N): newick}, ou
      - list [newick, newick, ...]
    Retour: {"D": np.ndarray, "labels": [...]}
    """
    # --- normalisation d'entrée + labels ---
    if isinstance(tree_data, dict):
        items  = sorted(tree_data.items(), key=lambda kv: (kv[0][0], kv[0][1]))
        keys   = [key for (key, _) in items]
        nwicks = [nw  for (_, nw) in items]
        if label_style in ("auto", "tuple"):
            labels = keys
        elif label_style == "string":
            labels = [f"cluster{k}-tree{n}" for (k, n) in keys]
        elif label_style == "index":
            labels = [f"T{i+1}" for i in range(len(keys))]
        elif label_style == "none":
            labels = []
        else:
            labels = keys
    else:
        # liste de newick : pas de (k,N)
        nwicks = list(tree_data)
        if label_style in ("auto", "index"):
            labels = [f"T{i+1}" for i in range(len(nwicks))]
        elif label_style == "none":
            labels = []
        else:
            labels = [f"T{i+1}" for i in range(len(nwicks))]

    n = len(nwicks)
    if n == 0:
        return {"D": np.zeros((0,0)), "labels": labels}

    # Matrices normalisées pour Weight/BranchLength/Height/Degree
    seqs = _collect_all_seqs(nwicks)
    Wn, BLn, Hn, DGn = _build_and_normalize(nwicks, seqs)

    # Pré-calc sets de feuilles et splits topologiques (speed-up)
    nodes_list, splits_list = _precompute_per_tree(nwicks)

    # Assemblage
    D = np.zeros((n, n), dtype=float)
    lam1, lam2, lam3, lam4, lam5 = lambdas

    for i in range(n):
        if progress and i % 8 == 0:
            print(f"[WMFD] {i}/{n}", flush=True)
        s1 = nodes_list[i]
        S1 = splits_list[i]
        for j in range(i+1, n):
            s2 = nodes_list[j]
            S2 = splits_list[j]
            Ic, Iu, C, U = _common_uncommon_indices(seqs, s1, s2)

            bl_c = _sum_abs_diff_rows(BLn, i, j, Ic, True)
            bl_u = _sum_abs_diff_rows(BLn, i, j, Iu, True)
            w_c  = _sum_abs_diff_rows(Wn,  i, j, Ic, True)
            w_u  = _sum_abs_diff_rows(Wn,  i, j, Iu, True)
            d_c  = _sum_abs_diff_rows(DGn, i, j, Ic, True)
            d_u  = _sum_abs_diff_rows(DGn, i, j, Iu, True)
            h_c  = _sum_abs_diff_rows(Hn,  i, j, Ic, True)
            h_u  = _sum_abs_diff_rows(Hn,  i, j, Iu, True)

            P  = _penalty(C, U)
            Usp = S1 | S2
            HD = 0.0 if len(Usp) == 0 else (1.0 - (len(S1 & S2) / len(Usp)))

            common   = lam1*bl_c + lam2*w_c + lam3*d_c + lam4*h_c
            uncommon = lam1*bl_u + lam2*w_u + lam3*d_u + lam4*h_u
            D[i, j] = D[j, i] = common + P*uncommon + lam5*HD

    np.fill_diagonal(D, 0.0)
    if progress: print(f"[WMFD] {n}/{n}", flush=True)
    return {"D": D, "labels": labels}

# --- test rapide ---
if __name__ == "__main__":

    newick_trees = [
        "((L1:0.6,L2:0.6):0.3,(L3:0.5,L4:0.5):0.2);",
        "(L1:0.7,(L2:0.5,L5:0.4):0.3);",
        "((L3:0.55,L4:0.55):0.25,L6:0.8);"
    ]
    out2 = compute_wmfd_all(newick_trees, progress=True, label_style="index")
    print("\n[LIST] labels:", out2["labels"])
    print("[LIST] D shape:", out2["D"].shape)
    print(out2["D"])
