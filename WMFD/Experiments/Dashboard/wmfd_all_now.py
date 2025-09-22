# wmfd_all_now.py  —  pipeline complet: génération -> WMFD -> clustering -> dashboard
# Dépendances: ete3, numpy, scikit-learn, matplotlib, pandas (facultatif)
from __future__ import annotations
import math, re, os, csv, random
from typing import List, Dict, Tuple, Set, Optional
import numpy as np
import matplotlib.pyplot as plt

import gptree_generate_structures as gen
from ete3 import Tree

# =============== PARAMS =================
MAX_RUNS: Optional[int] = None          # None = tous, sinon tronque pour debug (ex: 50)
MAX_TREES_PER_RUN = 128                 # échantillon par run pour accélérer
WRITE_FILES = False                     # True => export CSV par run (matrices)
OUT_DIR = "wmfd_out_all"

# Clustering/évaluation
AUTO_K = False          # False = k=K_true (recommandé pour benchmark WMFD)
RANGE_AROUND_K = 1      # si AUTO_K=True, testera {K-1, K, K+1} (borné 2..n-1)
RANDOM_STATE = 42

# =============== WMFD ===================

def _uniq_names(t: Tree) -> None:
    k = 0
    for lf in t.iter_leaves():
        if not lf.name:
            k += 1
            lf.name = f"__anon{k}"

def _precomp(t: Tree):
    _uniq_names(t)
    root = t.get_tree_root()
    BL: Dict[str, float] = {}
    H:  Dict[str, float] = {}
    W:  Dict[str, float] = {}
    D:  Dict[str, float] = {}
    Ls: Set[str] = set()
    Splits: Set[frozenset] = set()

    for lf in t.iter_leaves():
        nm = str(lf.name)
        Ls.add(nm)
        BL[nm] = float(lf.dist or 0.0)
        H[nm]  = float(root.get_distance(lf) or 0.0)
        W[nm]  = 1.0
        D[nm]  = float(len(lf.up.children) if lf.up else 1)

    for nd in t.traverse("postorder"):
        if not nd.is_leaf():
            cl = frozenset(le.name for le in nd.iter_leaves())
            if 1 < len(cl) < len(Ls):
                Splits.add(cl)

    return (t, BL, H, W, D, Ls, Splits)

def _nn(v1: float, v2: float, mn: float, mx: float) -> Tuple[float,float]:
    if mx > mn:
        return ((v1 - mn)/(mx - mn), (v2 - mn)/(mx - mn))
    return (0.0, 0.0)

def wmfd_pair(A, B, L1=0.30, L2=0.20, L3=0.25, L4=0.15, L5=0.10) -> float:
    (_t1, bl1, h1, w1, d1, Ls1, S1) = A
    (_t2, bl2, h2, w2, d2, Ls2, S2) = B

    U = Ls1 | Ls2
    I = Ls1 & Ls2
    TN = len(U); CN = len(I)

    blmin = min([bl1.get(x,0.0) for x in U] + [bl2.get(x,0.0) for x in U], default=0.0)
    blmax = max([bl1.get(x,0.0) for x in U] + [bl2.get(x,0.0) for x in U], default=0.0)
    hmin  = min([h1.get(x,0.0)  for x in U] + [h2.get(x,0.0)  for x in U], default=0.0)
    hmax  = max([h1.get(x,0.0)  for x in U] + [h2.get(x,0.0)  for x in U], default=0.0)
    dmin  = min([d1.get(x,0.0)  for x in U] + [d2.get(x,0.0)  for x in U], default=0.0)
    dmax  = max([d1.get(x,0.0)  for x in U] + [d2.get(x,0.0)  for x in U], default=0.0)

    cBL=cH=cW=cD = 0.0
    uBL=uH=uW=uD = 0.0
    for x in U:
        a,b = _nn(bl1.get(x,0.0), bl2.get(x,0.0), blmin, blmax); dBL = abs(a-b)
        a,b = _nn(h1.get(x,0.0),  h2.get(x,0.0),  hmin,  hmax);  dH  = abs(a-b)
        dW  = 0.0  # W constant
        a,b = _nn(d1.get(x,0.0),  d2.get(x,0.0),  dmin,  dmax);  dD  = abs(a-b)
        if x in I:
            cBL += dBL; cH += dH; cW += dW; cD += dD
        else:
            uBL += dBL; uH += dH; uW += dW; uD += dD

    CNs  = max(len(I), 1)
    UNCs = max(TN - len(I), 1)
    Wc = 0.30*(cBL/CNs) + 0.20*(cH/CNs) + 0.25*(cW/CNs) + 0.15*(cD/CNs)
    Wu = 0.30*(uBL/UNCs) + 0.20*(uH/UNCs) + 0.25*(uW/UNCs) + 0.15*(uD/UNCs)

    Usp = S1 | S2; Isp = S1 & S2
    HD  = 1.0 - (len(Isp)/len(Usp)) if len(Usp) > 0 else 0.0
    P   = 1.0 - (len(I) / TN) if TN > 0 else 0.0
    return float(P*Wu + Wc + 0.10*HD)

def wmfd_matrix_ete(trees: List[Tree]):
    feats = [_precomp(t) for t in trees]
    n = len(feats)
    D = np.zeros((n, n), dtype=float)
    for i in range(n):
        if i % 8 == 0:
            print(f"[all]   progress {i}/{n}", flush=True)
        for j in range(i+1, n):
            D[i, j] = D[j, i] = wmfd_pair(feats[i], feats[j])
    print(f"[all]   progress {n}/{n}", flush=True)
    return D

def _write_csv(stem: str, D: np.ndarray):
    os.makedirs(OUT_DIR, exist_ok=True)
    mpath = os.path.join(OUT_DIR, f"{stem}_wmfd_matrix.csv")
    with open(mpath, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        for i in range(D.shape[0]):
            w.writerow([f"{D[i,j]:.8f}" for j in range(D.shape[1])])
    print(f"[all]   CSV -> {mpath}", flush=True)

# =============== Génération via gptree (corrigé) ===============

print("[all] import générateur OK", flush=True)
runs = gen.generate_runs(
    Ks=[1,2,3,4], 
    Ls=[10,20,30,40,50,60,70],
    ns=[8,16],
    plevels=[0.3,0.5,0.7],
    noises=[0,25,50,75],          # maintenant UTILISÉ par le gptree
    reps=[0,1,2],
    return_format="ete",
    
)
print(f"[all] runs générés: {len(runs)}", flush=True)
if MAX_RUNS is not None:
    runs = runs[:MAX_RUNS]
    print(f"[all] runs limités à: {len(runs)}", flush=True)

all_mats: Dict[str, np.ndarray] = {}
for idx, r in enumerate(runs):
    trees = list(r.trees_ete or [])
    if not trees:
        print(f"[all] skip {r.run_name} (0 arbres)")
        continue
    if len(trees) > MAX_TREES_PER_RUN:
        trees = trees[:MAX_TREES_PER_RUN]
    print(f"[all] run {idx+1}/{len(runs)}: {r.run_name} | arbres utilisés={len(trees)}", flush=True)
    D = wmfd_matrix_ete(trees)
    all_mats[r.run_name] = D

    # stats rapides
    vals = D[np.triu_indices(D.shape[0], k=1)]
    mean = float(np.mean(vals)) if vals.size else 0.0
    print(f"[all]   mean pairwise WMFD = {mean:.6f}", flush=True)

    if WRITE_FILES:
        _write_csv(r.run_name, D)

print(f"[all] done ✅ | matrices en mémoire: {len(all_mats)}", flush=True)

