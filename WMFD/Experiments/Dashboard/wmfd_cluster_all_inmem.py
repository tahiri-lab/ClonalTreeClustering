#!/usr/bin/env python3
# -- coding: utf-8 --
from __future__ import annotations

import time
from typing import Dict, Set, Tuple, List, Optional
import numpy as np

from ete3 import Tree
from sklearn.metrics import adjusted_rand_score, silhouette_score
from sklearn.manifold import MDS
from sklearn.metrics import calinski_harabasz_score

print("[BOOT] importing generator…", flush=True)
import gptree_generate_structures as gen
print("[BOOT] generator imported ✓", flush=True)

# ===================== PARAMÈTRES GLOBAUX =====================
RANDOM_STATE = 42
AUTO_K      = True           # True = choisit K automatiquement (silhouette/CH)
CRITERION   = "silhouette"   # "silhouette" | "ch" | "both"
KMIN, KMAX  = 2, 10          # bornes pour la recherche de K si AUTO_K=True
MAX_TREES_PER_RUN: Optional[int] = None  # ex. 128 pour aller + vite (None = tout)

# Grilles "raisonnables" (à ajuster rapidement si besoin)
Ks      = [2, 3, 4]
Ls      = [20, 40, 60]
ns      = [16, 32]
plevels = [0.30, 0.50, 0.70]
noises  = [0, 25, 50]
reps    = [0]

# ===================== WMFD =====================
def _uniq_names(t: Tree) -> None:
    k = 0
    for lf in t.iter_leaves():
        if not lf.name:
            k += 1
            lf.name = f"__anon{k}"

def wmfd_precompute_tree(t: Tree):
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
        dW  = 0.0  # W constant (=1)
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

# ===================== K-MEDOIDS (dynMSC si dispo, sinon PAM) =====================
def _check_D(D: np.ndarray) -> np.ndarray:
    D = np.asarray(D, float)
    D = 0.5*(D + D.T)
    np.fill_diagonal(D, 0.0)
    D[D < 0] = 0.0
    return D

def pam_kmedoids(D: np.ndarray, k: int, max_iter: int = 60, seed: int = 42):
    rng = np.random.default_rng(seed)
    D = _check_D(D); n = D.shape[0]
    # BUILD
    costs = np.sum(D, axis=1)
    medoids = [int(np.argmin(costs))]
    while len(medoids) < k:
        cur = np.min(D[:, medoids], axis=1)
        best_gain, best_j = None, None
        for j in range(n):
            if j in medoids: continue
            new_near = np.minimum(cur, D[:, j])
            gain = float(np.sum(cur) - np.sum(new_near))
            if (best_gain is None) or (gain > best_gain + 1e-12):
                best_gain, best_j = gain, j
        if best_j is None:
            cand = [x for x in range(n) if x not in medoids]
            best_j = int(rng.choice(cand))
        medoids.append(int(best_j))
    # SWAP
    for _ in range(max_iter):
        improved = False
        current = float(np.sum(np.min(D[:, medoids], axis=1)))
        non = [i for i in range(n) if i not in medoids]
        for mi, m in enumerate(medoids):
            for h in non:
                trial = medoids.copy(); trial[mi] = h
                cost = float(np.sum(np.min(D[:, trial], axis=1)))
                if cost + 1e-12 < current:
                    medoids = trial; improved = True; break
            if improved: break
        if not improved: break
    labels = np.argmin(D[:, medoids], axis=1).astype(int)
    return np.array(medoids, int), labels

def choose_k_and_labels(D: np.ndarray, kmin=2, kmax=10, criterion="silhouette", seed=42):
    D = _check_D(D)

    # 1) dynMSC (si dispo) — uniquement pertinent pour "silhouette"
    if criterion in ("silhouette", "both"):
        try:
            import kmedoids  # pip install kmedoids
            dm = kmedoids.dynmsc(D, kmax, kmin)
            bestk = int(dm.bestk)
            sol = kmedoids.fasterpam(D, bestk, init='build')
            labels = np.asarray(sol.labels, int)
            medoids = np.asarray(sol.medoids, int)
            try:
                rangek = list(dm.rangek)
                sils   = list(dm.losses)  # silhouette par K
                sil_star = float(sils[rangek.index(bestk)]) if bestk in rangek else float("nan")
            except Exception:
                rangek, sils, sil_star = list(range(kmin, kmax+1)), [], float("nan")

            if criterion == "silhouette":
                return bestk, labels, medoids, rangek, sils, sil_star, None

            # if criterion == "both", on calcule aussi CH via MDS
            X = MDS(n_components=4, dissimilarity="precomputed", random_state=seed,
                    n_init=1, max_iter=200).fit_transform(D)
            ch_star = float(calinski_harabasz_score(X, labels))
            # pour garder la trace des CH par K on refait un petit tour
            chs = []
            for k in rangek:
                if k < 2 or k >= D.shape[0]:
                    chs.append(np.nan)
                    continue
                _, lab_k = pam_kmedoids(D, k=k, seed=seed)
                try:
                    chs.append(float(calinski_harabasz_score(X, lab_k)))
                except Exception:
                    chs.append(np.nan)
            return bestk, labels, medoids, rangek, sils, sil_star, chs

        except Exception as e:
            print(f"[dynMSC] indisponible ({e}). Fallback PAM+{criterion}.", flush=True)

    # 2) Fallback full-Python (PAM + métrique choisie)
    best_k, best_key, best_labels = None, -1e9, None
    rangek = list(range(max(2, kmin), max(3, kmax)+1))
    sils: List[float] = []
    chs : List[float] = []

    # Embedding MDS si CH utilisé
    X = None
    if criterion in ("ch", "both"):
        X = MDS(n_components=4, dissimilarity="precomputed",
                random_state=seed, n_init=1, max_iter=200).fit_transform(D)

    for k in rangek:
        if k < 2 or k >= D.shape[0]: 
            if criterion in ("ch", "both"): chs.append(np.nan)
            if criterion in ("silhouette", "both"): sils.append(np.nan)
            continue
        _, labels = pam_kmedoids(D, k=k, seed=seed)
        key = 0.0
        s_val = np.nan; c_val = np.nan
        if criterion in ("silhouette", "both"):
            try:
                s_val = float(silhouette_score(D, labels, metric="precomputed"))
            except Exception:
                s_val = float("-inf")
            sils.append(s_val)
            key += s_val
        if criterion in ("ch", "both"):
            try:
                c_val = float(calinski_harabasz_score(X, labels))
            except Exception:
                c_val = float("-inf")
            chs.append(c_val)
            # faible poids ajouté pour départager si égalité de silhouette
            if criterion == "both":
                key += 1e-3 * (c_val if np.isfinite(c_val) else 0.0)

        if key > best_key:
            best_k, best_key, best_labels = k, key, labels

    medoids, best_labels = pam_kmedoids(D, k=best_k, seed=seed)
    sil_star = (np.nanmax(sils) if sils else None) if criterion in ("silhouette","both") else None
    return best_k, best_labels, medoids, rangek, sils if sils else None, sil_star, chs if chs else None

# ===================== PIPELINE TOUT-EN-UN =====================
def run_all():
    print("[ALL] génération de tous les runs…", flush=True)
    runs = gen.generate_runs(
        Ks=Ks, Ls=Ls, ns=ns, plevels=plevels, noises=noises,
        reps=reps, return_format="ete"
    )
    print(f"[ALL] runs générés: {len(runs)}", flush=True)

    results: List[Dict[str, object]] = []
    t0 = time.perf_counter()

    for idx, r in enumerate(runs, 1):
        trees = list(r.trees_ete or [])
        if not trees:
            print(f"[RUN {idx}] skip {r.run_name} (0 arbres)", flush=True)
            continue
        if MAX_TREES_PER_RUN is not None and len(trees) > MAX_TREES_PER_RUN:
            trees = trees[:MAX_TREES_PER_RUN]

        print(f"\n[RUN {idx}/{len(runs)}] {r.run_name} | arbres utilisés={len(trees)} | K_true={r.K}", flush=True)

        # (1) WMFD
        D = wmfd_matrix_ete(trees, progress=True)
        vals = D[np.triu_indices(D.shape[0], k=1)]
        print(f"[WMFD] shape={D.shape} | min={vals.min():.4f} | max={vals.max():.4f} | mean={vals.mean():.4f}", flush=True)

        # (2) Clustering
        if AUTO_K:
            k_star, y_pred, medoids, rangek, sils, sil_star, chs = choose_k_and_labels(
                D, kmin=KMIN, kmax=min(KMAX, D.shape[0]-1), criterion=CRITERION, seed=RANDOM_STATE
            )
        else:
            k_star = r.K
            medoids, y_pred = pam_kmedoids(D, k=r.K, seed=RANDOM_STATE)
            rangek, sils, chs, sil_star = [], None, None, None

        # (3) ARI + PRINTS CLAIRS —
        from collections import defaultdict
        # --- y_pred_full : identifiant cluster_index pour chaque arbre prédit ---
        y_pred = np.asarray(y_pred, dtype=int)  # au cas où
        ctr = defaultdict(int)
        y_pred_full = []
        for lbl in y_pred:
         c = int(lbl)
         idx = ctr[lbl]
        y_pred_full.append(f"{lbl}_{idx}")
        ctr[lbl]+=1

        y_true_full = list(r.true_labels)[:len(y_pred)]  
        def _to_cluster_id(s):
            try:
                return int(str(s).split("_", 1)[0])
            except Exception:
                try:
                    return int(s)
                except Exception:
                    return -1
        y_true_cluster = np.array([_to_cluster_id(s) for s in y_true_full], dtype=int)

        # ARI (clusters only)
        ari = float(adjusted_rand_score(y_true_cluster, y_pred))

        # tailles par cluster (vrai / prédit)
        from collections import Counter
        def _sizes(lbls):
            c = Counter(lbls)
            return "{" + ", ".join([f"{k}:{c[k]}" for k in sorted(c)]) + "}"
        sizes_true = _sizes(y_true_cluster.tolist())
        sizes_pred = _sizes(y_pred.tolist())

        # silhouette* (si déjà calculée), sinon recalcul rapide
        if sils is not None and len(sils) > 0 and np.isfinite(np.nanmax(sils)):
            sil_star_val = float(sil_star)
        else:
            try:
                sil_star_val = float(silhouette_score(D, y_pred, metric="precomputed"))
            except Exception:
                sil_star_val = float("nan")

        # médoines (affichage sûr) + leur cluster vrai
        medoids_list = [] if medoids is None else [int(x) for x in np.asarray(medoids).ravel().tolist()]
        med_true = [int(y_true_cluster[m]) for m in medoids_list] if medoids_list else []

        # (optionnel) Purity
        def purity(y_true_int, y_pred_int):
            from collections import defaultdict
            groups = defaultdict(list)
            for yt, yp in zip(y_true_int, y_pred_int):
                groups[yp].append(yt)
            return sum(max(Counter(g).values()) for g in groups.values()) / len(y_true_int)
        pur = purity(y_true_cluster, y_pred)

        # DEBUG prints clairs
        print(f"[DEBUG] y_true_full[:20]={y_true_full[:20]}")
        print(f"[DEBUG] sizes_true={sizes_true} | sizes_pred={sizes_pred}")
        print(f"[DEBUG] y_pred_full[:20]={y_pred_full[:20]}")
        if medoids_list:
            print(f"[DEBUG] medoids={medoids_list} | medoids_true_clusters={med_true}")
        print(f"[METRICS] ARI_cluster={ari:.3f} | Silhouette*={sil_star_val:.3f} | Purity={pur:.3f}")

        # (4) Stocke pour résumé
        results.append({
            "name": r.run_name,
            "K_true": r.K,
            "L": r.L,
            "n_per_group": r.n_per_group,
            "plevel": r.plevel,
            "noise": r.noise_pct,
            "N": len(trees),
            "k_star": int(k_star),
            "ari": float(ari),
            "sil_star": float(sil_star_val) if np.isfinite(sil_star_val) else np.nan,
            "purity": float(pur),
        })

    # --------- Résumé lisible : ARI moyen par niveau de bruit ---------
    if results:
        print("\n[SUMMARY] ARI moyen par noise (%):", flush=True)
        noises_sorted = sorted({r["noise"] for r in results})
        for p in noises_sorted:
            aris = [r["ari"] for r in results if r["noise"] == p]
            if aris:
                print(f"  noise={int(p):2d}% -> mean ARI={np.mean(aris):.3f}", flush=True)

        # Optionnel : top/bas runs
        res_sorted = sorted(results, key=lambda x: x["ari"], reverse=True)
        print("\n[TOP 5] meilleurs ARI:", flush=True)
        for row in res_sorted[:5]:
            print(f"  {row['name']} | K*={row['k_star']} | ARI={row['ari']:.3f}", flush=True)
        print("\n[BOTTOM 5] pires ARI:", flush=True)
        for row in res_sorted[-5:]:
            print(f"  {row['name']} | K*={row['k_star']} | ARI={row['ari']:.3f}", flush=True)

# ===================== MAIN =====================
if __name__ == "__main__":
    run_all()