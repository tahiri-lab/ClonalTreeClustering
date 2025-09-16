# -*- coding: utf-8 -*-
"""
Step 1c) Clustering from pairwise distances (WMFD or RF) -> predicted partitions.

Entry:
- -matrix: CSV N N in the format:
    id,1.1,1.2,1.3,...
    1.1,0,  d12, d13,...
    1.2,d21, 0,  d23,...
    ...
  (compatible with wmfd_matrix.csv generated in 1b)

Optional:
- --labels: labels.csv (id,true_cluster) to calculate ARI (1d as a bonus)

Outputs:
- clusters_pred.csv   : id,pred_cluster   (predicted partition)
- medoids.csv   : cluster,medoid_id
- dist_to_medoids.csv: id,<medoid1>,...,<medoidK>  (distance from each tree to each medoid)
- metrics.csv   : k,objective,sum_intra,ARI(optional),n_iter,n_init,seed

Algo: K-medoids (simplified PAM) on pre-calculated distances, with n_init random restarts.
"""

import csv
import argparse
import random
from typing import List, Tuple

# --------- Read/Write Utilities ----------

def read_distance_matrix(path: str):
    with open(path, newline="", encoding="utf-8") as f:
        rdr = csv.reader(f)
        rows = list(rdr)
    header = rows[0]
    ids = header[1:]
    n = len(ids)
    D = [[0.0]*n for _ in range(n)]
    for i, row in enumerate(rows[1:]):
        rid = row[0]
        if rid != ids[i]:
            # we tolerate but we force the alignment by position
            pass
        vals = row[1:]
        if len(vals) != n:
            raise ValueError(f"Ligne {i+1}: taille attendue {n}, trouvée {len(vals)}")
        for j, v in enumerate(vals):
            D[i][j] = float(v)
    # balance + diag=0
    for i in range(n):
        D[i][i] = 0.0
        for j in range(i+1, n):
            m = 0.5*(D[i][j] + D[j][i])
            D[i][j] = D[j][i] = m
    return ids, D

def write_clusters(path: str, ids: List[str], labels: List[int]):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["id","pred_cluster"])
        for i, lab in zip(ids, labels):
            w.writerow([i, lab])

def write_medoids(path: str, medoids_idx: List[int], ids: List[str]):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["cluster","medoid_id"])
        for c, m in enumerate(medoids_idx, start=1):
            w.writerow([c, ids[m]])

def write_dist_to_medoids(path: str, ids: List[str], medoids_idx: List[int], D):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        header = ["id"] + [f"medoid_{c+1}" for c in range(len(medoids_idx))]
        w.writerow(header)
        for i, iid in enumerate(ids):
            row = [iid] + [f"{D[i][m]:.6f}" for m in medoids_idx]
            w.writerow(row)

def write_metrics(path: str, k: int, obj: float, n_iter: int, n_init: int, seed: int, ari: float = None):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["k","objective","n_iter","n_init","seed","ARI"])
        w.writerow([k, f"{obj:.6f}", n_iter, n_init, seed, (f"{ari:.6f}" if ari is not None else "")])

# --------- ARI (optional for 1d, here as a bonus if labels provided) ---------

def _comb2(x): return x*(x-1)/2.0

def adjusted_rand_index(y_true: List[int], y_pred: List[int]) -> float:
    from collections import defaultdict
    n = len(y_true)
    by_t = defaultdict(list); by_p = defaultdict(list)
    for i,(t,p) in enumerate(zip(y_true, y_pred)):
        by_t[t].append(i); by_p[p].append(i)
    nij_sum = 0.0
    for ti in by_t.values():
        set_ti = set(ti)
        for pj in by_p.values():
            nij = len(set_ti & set(pj))
            nij_sum += _comb2(nij)
    ai_sum = sum(_comb2(len(v)) for v in by_t.values())
    bj_sum = sum(_comb2(len(v)) for v in by_p.values())
    n_pairs = _comb2(n)
    expected = (ai_sum * bj_sum) / n_pairs if n_pairs else 0.0
    max_index = 0.5*(ai_sum + bj_sum)
    denom = max_index - expected
    if denom == 0:
        return 1.0 if nij_sum == max_index else 0.0
    return (nij_sum - expected) / denom

def read_true_labels(path: str, ids: List[str]):
    # labels.csv: id,true_cluster
    m = {}
    with open(path, newline="", encoding="utf-8") as f:
        r = csv.DictReader(f)
        for row in r:
            m[row["id"]] = int(row["true_cluster"])
    return [m[i] for i in ids]

# ---------- K-medoids ----------

def assign_labels(D, medoids):
    n = len(D)
    labels = [0]*n
    for i in range(n):
        best, best_c = float("inf"), None
        for c, m in enumerate(medoids):
            d = D[i][m]
            if d < best:
                best, best_c = d, c+1  # clusters 1..k
        labels[i] = best_c
    return labels

def update_medoids(D, labels, k):
    n = len(D)
    medoids = [None]*k
    for c in range(1, k+1):
        members = [i for i in range(n) if labels[i] == c]
        if not members:
            continue
        best_sum, best_i = float("inf"), members[0]
        for i in members:
            s = 0.0
            for j in members:
                s += D[i][j]
            if s < best_sum:
                best_sum, best_i = s, i
        medoids[c-1] = best_i
    return medoids

def objective(D, labels, medoids):
    # sum of the distances to the medoid of each cluster
    s = 0.0
    for i, lab in enumerate(labels):
        m = medoids[lab-1]
        s += D[i][m]
    return s

def k_medoids(D, k, seed=0, n_init=10, max_iter=100):
    rng = random.Random(seed)
    n = len(D)
    best_obj, best_labels, best_meds, best_iters = float("inf"), None, None, 0
    for run in range(n_init):
        # random init without duplicate
        medoids = rng.sample(range(n), k)
        labels = assign_labels(D, medoids)
        it = 0
        while it < max_iter:
            it += 1
            new_medoids = update_medoids(D, labels, k)
            # if medoid empty, we resample
            for idx, m in enumerate(new_medoids):
                if m is None:
                    cand = rng.choice([i for i in range(n) if i not in new_medoids and i not in medoids])
                    new_medoids[idx] = cand
            if set(new_medoids) == set(medoids):
                break
            medoids = new_medoids
            labels = assign_labels(D, medoids)
        obj = objective(D, labels, medoids)
        if obj < best_obj:
            best_obj, best_labels, best_meds, best_iters = obj, labels[:], medoids[:], it
    return best_labels, best_meds, best_obj, best_iters

# ---------- Main ----------

def main():
    ap = argparse.ArgumentParser(description="1c) Clustering (K-medoids) from pairwise distance matrix")
    ap.add_argument("--matrix", type=str, required=True, help="wmfd_matrix.csv (ou RF matrix)")
    ap.add_argument("--k", type=int, required=True, help="nombre de clusters")
    ap.add_argument("--labels", type=str, default=None, help="labels.csv pour ARI (optionnel)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--n_init", type=int, default=20, help="redémarrages aléatoires")
    ap.add_argument("--max_iter", type=int, default=100)
    ap.add_argument("--out_pred", type=str, default="clusters_pred.csv")
    ap.add_argument("--out_medoids", type=str, default="medoids.csv")
    ap.add_argument("--out_kxn", type=str, default="dist_to_medoids.csv")
    ap.add_argument("--out_metrics", type=str, default="metrics.csv")
    args = ap.parse_args()

    ids, D = read_distance_matrix(args.matrix)

    labels, medoids, obj, n_iter = k_medoids(D, args.k, seed=args.seed, n_init=args.n_init, max_iter=args.max_iter)

    # outputs
    write_clusters(args.out_pred, ids, labels)
    write_medoids(args.out_medoids, medoids, ids)
    write_dist_to_medoids(args.out_kxn, ids, medoids, D)

    ari = None
    if args.labels:
        y_true = read_true_labels(args.labels, ids)
        ari = adjusted_rand_index(y_true, labels)

    write_metrics(args.out_metrics, args.k, obj, n_iter, args.n_init, args.seed, ari)
    print("Done.")
    if ari is not None:
        print(f"ARI = {ari:.6f}")
    print(f"Objective(sum intra) = {obj:.6f}")

if __name__ == "__main__":
    main()
