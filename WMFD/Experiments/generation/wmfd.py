#!/usr/bin/env python3
# -- coding: utf-8 --
from __future__ import annotations

import time
from typing import Dict, Set, Tuple, List, Optional
import numpy as np

from ete3 import Tree

print("[BOOT] importing generator…", flush=True)
import gptree_generate_structures as gen
print("[BOOT] generator imported ✓", flush=True)

# ===================== PARAMÈTRES GLOBAUX =====================

DEFAULT_Ks      = (1, 2, 3, 4)
DEFAULT_Ls      = (10, 20, 30, 40, 50, 60, 70)
DEFAULT_ns      = (8 , 16)
DEFAULT_plevels = (0.30, 0.50, 0.70)
DEFAULT_noises  = (0, 25, 50, 75 )
DEFAULT_reps    = (0,1)

# ===================== WMFD =====================
def _uniq_names(t: Tree) -> None:
    k = 0
    for lf in t.iter_leaves():
        if not lf.name:
            k += 1
            lf.name = f"__anon{k}"

def wmfd_precompute_tree(t: Tree):
    """
    Pré-calcul pour WMFD : longueurs de branches terminales, hauteurs, etc.
    Retourne un tuple (t, BL, H, W, D, leaf_set, splits).
    """
    _uniq_names(t)
    root = t.get_tree_root()
    BL: Dict[str, float] = {}
    H:  Dict[str, float] = {}
    W:  Dict[str, float] = {}
    D:  Dict[str, float] = {}
    for lf in t.iter_leaves():
        nm = str(lf.name)
        BL[nm] = float(lf.dist or 0.0)
        H[nm]  = float(root.get_distance(lf) or 0.0)
        W[nm]  = 1.0
        D[nm]  = float(len(lf.up.children) if lf.up else 1)
    leaf_set: Set[str] = set(BL.keys())
    splits: Set[frozenset] = set()
    for nd in t.traverse("postorder"):
        if nd.is_leaf():
            continue
        clade = frozenset(le.name for le in nd.iter_leaves())
        if 1 < len(clade) < len(leaf_set):
            splits.add(clade)
    return (t, BL, H, W, D, leaf_set, splits)

def _pair_norm(v1: float, v2: float, mn: float, mx: float) -> Tuple[float, float]:
    if mx > mn:
        return ((v1 - mn) / (mx - mn), (v2 - mn) / (mx - mn))
    return (0.0, 0.0)

def wmfd_pair(A, B, L1=0.30, L2=0.20, L3=0.25, L4=0.15, L5=0.10) -> float:
    """
    Distance WMFD entre deux arbres pré-calculés A et B.
    """
    (_t1, bl1, h1, w1, d1, Ls1, S1) = A
    (_t2, bl2, h2, w2, d2, Ls2, S2) = B
    U = Ls1 | Ls2
    I = Ls1 & Ls2
    TN = len(U); CN = len(I)

    bl_min = min([bl1.get(x,0.0) for x in U] + [bl2.get(x,0.0) for x in U], default=0.0)
    bl_max = max([bl1.get(x,0.0) for x in U] + [bl2.get(x,0.0) for x in U], default=0.0)
    h_min  = min([h1.get(x,0.0)  for x in U] + [h2.get(x,0.0)  for x in U], default=0.0)
    h_max  = max([h1.get(x,0.0)  for x in U] + [h2.get(x,0.0)  for x in U], default=0.0)
    d_min  = min([d1.get(x,0.0)  for x in U] + [d2.get(x,0.0)  for x in U], default=0.0)
    d_max  = max([d1.get(x,0.0)  for x in U] + [d2.get(x,0.0)  for x in U], default=0.0)

    cBL=cH=cW=cD = 0.0
    uBL=uH=uW=uD = 0.0
    for x in U:
        a,b = _pair_norm(bl1.get(x,0.0), bl2.get(x,0.0), bl_min, bl_max); dBL = abs(a-b)
        a,b = _pair_norm(h1.get(x,0.0),  h2.get(x,0.0),  h_min,  h_max);  dH  = abs(a-b)
        dW  = 0.0  
        a,b = _pair_norm(d1.get(x,0.0),  d2.get(x,0.0),  d_min,  d_max);  dD  = abs(a-b)
        if x in I:
            cBL += dBL; cH += dH; cW += dW; cD += dD
        else:
            uBL += dBL; uH += dH; uW += dW; uD += dD

    CNs  = max(CN, 1)
    UNCs = max(TN - CN, 1)
    Wc = L1*(cBL/CNs) + L2*(cH/CNs) + L3*(cW/CNs) + L4*(cD/CNs)
    Wu = L1*(uBL/UNCs) + L2*(uH/UNCs) + L3*(uW/UNCs) + L4*(uD/UNCs)

    Usp = S1 | S2; Isp = S1 & S2
    HD  = 1.0 - (len(Isp)/len(Usp)) if len(Usp) > 0 else 0.0
    P   = 1.0 - (CN / TN) if TN > 0 else 0.0
    return float(P*Wu + Wc + L5*HD)

def wmfd_matrix_ete(trees: List[Tree], progress: bool = True) -> np.ndarray:
    """Construit la matrice de distances WMFD (pairwise) sur une liste d’arbres ETE."""
    feats = [wmfd_precompute_tree(t) for t in trees]
    n = len(feats)
    D = np.zeros((n, n), dtype=float)
    for i in range(n):
        if progress and (i % 8 == 0):
            print(f"[WMFD] {i}/{n}", flush=True)
        fi = feats[i]
        for j in range(i+1, n):
            D[i, j] = D[j, i] = wmfd_pair(fi, feats[j])
    if progress:
        print(f"[WMFD] {n}/{n}", flush=True)
    np.fill_diagonal(D, 0.0)
    return D

# ===================== PIPELINE IN-MEMORY =====================
def compute_all_wmfd_inmem(
    Ks: Tuple[int, ...] = DEFAULT_Ks,
    Ls: Tuple[int, ...] = DEFAULT_Ls,
    ns: Tuple[int, ...] = DEFAULT_ns,
    plevels: Tuple[float, ...] = DEFAULT_plevels,
    noises: Tuple[float, ...] = DEFAULT_noises,
    reps: Tuple[int, ...] = DEFAULT_reps,
    return_format: str = "ete",
    progress: bool = True,
):
    """
      out[run_name] = {
          "D": np.ndarray (matrice de distances pairwise),
          "trees": List[ete3.Tree],
          "meta": {"K", "L", "n_per_group", "plevel", "noise", "rep"}
      }
    """
    print("[ALL] génération…", flush=True)
    runs = gen.generate_runs(Ks, Ls, ns, plevels, noises, reps, return_format=return_format)
    print(f"[ALL] runs générés: {len(runs)}", flush=True)

    out = {}
    t0 = time.perf_counter()
    for i, r in enumerate(runs, 1):
        trees = list(r.trees_ete or [])
        if not trees:
            print(f"[SKIP] {r.run_name} (0 arbres)", flush=True)
            continue

        if progress:
            print(f"\n[RUN {i}/{len(runs)}] {r.run_name} | N={len(trees)}", flush=True)

        D = wmfd_matrix_ete(trees, progress=progress)
        vals = D[np.triu_indices(D.shape[0], 1)]
        if progress:
            print(f"[WMFD] shape={D.shape} | min={vals.min():.4f} | max={vals.max():.4f} | mean={vals.mean():.4f}", flush=True)

        out[r.run_name] = {
            "D": D,
            "trees": trees,
            "meta": {
                "K": r.K,
                "L": r.L,
                "n_per_group": r.n_per_group,
                "plevel": r.plevel,
                "noise": r.noise_pct,
                "rep": r.rep,
            },
            "labels_true": list(r.true_labels)
        }
    t1 = time.perf_counter()
    print(f"\n[ALL] terminé ✓ | temps total: {t1 - t0:.1f}s | objets en mémoire: {len(out)}", flush=True)
    return out

# ===================== MAIN =====================
if __name__ == "__main__":
    data = compute_all_wmfd_inmem()

    if data:
        first_run = next(iter(data))
        D0 = data[first_run]["D"]
        print(f"[CHECK] 1er run: {first_run} | D shape={D0.shape}", flush=True)