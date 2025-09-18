# step1c_all_runs.py
# -- coding: utf-8 --
"""
1c) K-Medoids from WMFD matrices — zero-argument, scans all runs automatically.

Behavior
--------
- Start roots: current working directory AND ./A_noise_dense (if exists)
- Discover runs by presence of a "wmfd" subfolder
- For each run, pick matrix:
    prefer clean: wmfd/trees_wmfd_matrix.csv
    fallback    : wmfd/trees_noisy_wmfd_matrix.csv
- Try K in 2..min(12, N-1) (K=1 skipped because silhouette is undefined)
- Select best K by:
    1) highest Silhouette (metric='precomputed')
    2) tie-break: highest Calinski–Harabasz (medoid adaptation)
    3) tie-break: lowest Davies–Bouldin (medoid adaptation)
- If labels are available (trees.csv or labels.csv), compute ARI (not used for selection)
- Write per-run outputs under: run_dir/cluster_kmedoids/K={best}/
- Write global summary: ./kmedoids_summary.csv

Notes
-----
- Uses sklearn-extra KMedoids if available; otherwise falls back to a pure-Python PAM.
- No CLI arguments: just run python step1c_all_runs.py.
"""

import os
import csv
import math
from typing import List, Tuple, Optional

import numpy as np
from sklearn.metrics import silhouette_score

# Try sklearn-extra
_BACKEND = None
try:
    from sklearn_extra.cluster import KMedoids as _SKEX_KMedoids  # pip install scikit-learn-extra
    _BACKEND = "sklearn-extra"
except Exception:
    _BACKEND = None


# ----------------------------- I/O helpers -----------------------------

def _to_float(x: str) -> float:
    x = (x or "").strip().replace(",", ".")
    return float(x) if x else 0.0

def read_distance_matrix(path: str) -> Tuple[List[str], np.ndarray]:
    with open(path, newline="", encoding="utf-8") as f:
        rdr = csv.reader(f)
        rows = list(rdr)
    if not rows:
        raise ValueError(f"Empty matrix file: {path}")
    header = [h.strip() for h in rows[0]]
    ids = header[1:]
    n = len(ids)
    D = np.zeros((n, n), dtype=float)
    for i, row in enumerate(rows[1:]):
        vals = row[1:]
        if len(vals) != n:
            raise ValueError(f"{os.path.basename(path)} row {i+1}: expected {n} values, got {len(vals)}")
        for j, v in enumerate(vals):
            D[i, j] = _to_float(v)
    D = 0.5 * (D + D.T)
    np.fill_diagonal(D, 0.0)
    return ids, D

def read_true_labels(path: str, ids: List[str]) -> Optional[List[int]]:
    if not os.path.exists(path):
        return None
    with open(path, newline="", encoding="utf-8") as f:
        rdr = csv.DictReader(f)
        if not rdr.fieldnames:
            return None
        keys = {k: (k or "").strip().lower() for k in rdr.fieldnames}
        id_key = next((k for k, lk in keys.items() if lk == "id"), None)
        y_key  = next((k for k, lk in keys.items() if lk in ("true_cluster","true","label","y","cluster")), None)
        if not id_key or not y_key:
            return None
        m = {}
        for row in rdr:
            try:
                m[row[id_key].strip()] = int(float(str(row[y_key]).strip().replace(",", ".")))
            except Exception:
                pass
    try:
        return [m[i] for i in ids]
    except KeyError:
        return None

def write_clusters(path: str, ids: List[str], labels_1based: List[int]) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f); w.writerow(["id", "pred_cluster"])
        for i, lab in zip(ids, labels_1based):
            w.writerow([i, lab])

def write_medoids(path: str, medoids_idx: List[int], ids: List[str]) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f); w.writerow(["cluster", "medoid_id"])
        for c, m in enumerate(medoids_idx, start=1):
            w.writerow([c, ids[m]])

def write_dist_to_medoids(path: str, ids: List[str], medoids_idx: List[int], D: np.ndarray) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        header = ["id"] + [f"medoid_{c+1}" for c in range(len(medoids_idx))]
        w.writerow(header)
        cols = D[:, medoids_idx]
        for i, iid in enumerate(ids):
            row = [iid] + [f"{d:.6f}" for d in cols[i, :]]
            w.writerow(row)

def write_metrics(path: str, row: dict) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    header = ["k","objective","n_iter","n_init","seed","ARI","silhouette","calinski_harabasz","davies_bouldin","backend"]
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=header)
        w.writeheader()
        w.writerow(row)


# ----------------------------- Label metrics -----------------------------

def adjusted_rand_index(y_true: List[int], y_pred: List[int]) -> float:
    from collections import defaultdict
    n = len(y_true)
    if n == 0: return float("nan")
    by_t = defaultdict(list); by_p = defaultdict(list)
    for idx, (t, p) in enumerate(zip(y_true, y_pred)):
        by_t[t].append(idx); by_p[p].append(idx)
    def comb2(x): return x * (x - 1) / 2.0
    nij_sum = 0.0
    for ti in by_t.values():
        set_ti = set(ti)
        for pj in by_p.values():
            nij = len(set_ti & set(pj))
            nij_sum += comb2(nij)
    ai_sum = sum(comb2(len(v)) for v in by_t.values())
    bj_sum = sum(comb2(len(v)) for v in by_p.values())
    n_pairs = comb2(n)
    expected = (ai_sum * bj_sum) / n_pairs if n_pairs else 0.0
    max_index = 0.5 * (ai_sum + bj_sum)
    denom = max_index - expected
    if denom == 0:
        return 1.0 if nij_sum == max_index else 0.0
    return (nij_sum - expected) / denom

def silhouette_from_precomputed(D: np.ndarray, labels_1: List[int]) -> float:
    labs0 = np.asarray(labels_1, dtype=int) - 1
    if len(np.unique(labs0)) < 2:
        return float("nan")
    return float(silhouette_score(D, labs0, metric="precomputed"))

def objective_sum_to_medoids(D: np.ndarray, labels_1: List[int], medoids: List[int]) -> float:
    idx = np.arange(D.shape[0])
    meds = np.array([medoids[lab - 1] for lab in labels_1], dtype=int)
    return float(np.sum(D[idx, meds]))

def ch_medoids(D: np.ndarray, labels_1: List[int], medoids: List[int]) -> float:
    n = D.shape[0]
    k = len(set(labels_1))
    if k < 2 or n <= k:
        return float("nan")
    within = 0.0
    labs = np.asarray(labels_1, dtype=int)
    for c in range(1, k + 1):
        members = np.where(labs == c)[0]
        if members.size == 0:
            return float("nan")
        m = medoids[c - 1]
        within += float(np.sum(D[members, m] ** 2))
    sums = np.sum(D, axis=1)
    g = int(np.argmin(sums))
    total = float(np.sum(D[:, g] ** 2))
    between = max(total - within, 0.0)
    return (between / (k - 1)) / (within / (n - k)) if within > 0 else float("inf")

def db_medoids(D: np.ndarray, labels_1: List[int], medoids: List[int]) -> float:
    k = len(set(labels_1))
    if k < 2:
        return float("nan")
    labs = np.asarray(labels_1, dtype=int)
    S = np.zeros(k, dtype=float)
    for c in range(1, k + 1):
        members = np.where(labs == c)[0]
        if members.size == 0:
            return float("nan")
        m = medoids[c - 1]
        S[c - 1] = float(np.mean(D[members, m])) if members.size > 0 else 0.0
    R = np.zeros((k, k), dtype=float)
    for i in range(k):
        for j in range(k):
            if i == j:
                R[i, j] = 0.0
            else:
                dij = float(D[medoids[i], medoids[j]])
                R[i, j] = (S[i] + S[j]) / dij if dij > 0 else np.inf
    Ri = np.max(R, axis=1)
    return float(np.mean(Ri))


# ----------------------------- K-Medoids engines -----------------------------

def _assign_labels(D: np.ndarray, medoids: List[int]) -> List[int]:
    d_to_meds = D[:, medoids]
    return list(np.argmin(d_to_meds, axis=1) + 1)

def _update_medoids(D: np.ndarray, labels_1: List[int], k: int) -> List[int]:
    medoids = [None] * k
    labs = np.asarray(labels_1, dtype=int)
    for c in range(1, k + 1):
        members = np.where(labs == c)[0]
        if members.size == 0:
            medoids[c - 1] = None
            continue
        sub = D[np.ix_(members, members)]
        sums = sub.sum(axis=1)
        medoids[c - 1] = int(members[int(np.argmin(sums))])
    return medoids

def pam_fallback(D: np.ndarray, k: int, seed=0, n_init=10, max_iter=100):
    rng = np.random.default_rng(seed)
    n = D.shape[0]
    best = (math.inf, None, None, 0)
    for r in range(max(1, n_init)):
        medoids = rng.choice(n, size=k, replace=False).tolist()
        labels = _assign_labels(D, medoids)
        it = 0
        while it < max_iter:
            it += 1
            new_meds = _update_medoids(D, labels, k)
            # repair empty clusters
            for i, m in enumerate(new_meds):
                if m is None:
                    dmin = np.min(D[:, medoids], axis=1)
                    cand_pool = [j for j in range(n) if j not in new_meds and j not in medoids]
                    cand = int(max(cand_pool, key=lambda j: dmin[j])) if cand_pool else int(rng.integers(0, n))
                    new_meds[i] = cand
            if set(new_meds) == set(medoids):
                break
            medoids = new_meds
            labels = _assign_labels(D, medoids)
        obj = objective_sum_to_medoids(D, labels, medoids)
        if obj < best[0]:
            best = (obj, labels[:], medoids[:], it)
    return best[1], best[2], best[0], best[3]

def run_kmedoids(D: np.ndarray, k: int, seed=0, n_init=10, max_iter=100):
    if k == 1:
        sums = D.sum(axis=1)
        m = int(np.argmin(sums))
        return [1]*D.shape[0], [m], float(np.sum(D[:, m])), 0, "k=1 closed-form"
    if _BACKEND == "sklearn-extra":
        best = (math.inf, None, None, 0, "")
        for r in range(max(1, n_init)):
            rs = seed + r
            km = _SKEX_KMedoids(
                n_clusters=k, metric="precomputed", method="pam",
                init="random", max_iter=max_iter, random_state=rs
            )
            km.fit(D)
            labels = list(km.labels_.astype(int) + 1)
            medoids = list(map(int, km.medoid_indices_))
            obj = objective_sum_to_medoids(D, labels, medoids)
            iters = getattr(km, "n_iter_", 0)
            if obj < best[0]:
                best = (obj, labels, medoids, iters, "sklearn-extra KMedoids(PAM)")
        return best[1], best[2], best[0], best[3], best[4]
    else:
        labels, medoids, obj, iters = pam_fallback(D, k, seed=seed, n_init=10, max_iter=max_iter)
        return labels, medoids, obj, iters, "fallback PAM (pure Python)"


# ----------------------------- K selection -----------------------------

def choose_best_K(D: np.ndarray, K_list: List[int], seed=0, max_iter=100):
    results = []
    for k in K_list:
        labels, medoids, obj, iters, backend = run_kmedoids(D, k, seed=seed, n_init=10, max_iter=max_iter)
        sil = silhouette_from_precomputed(D, labels) if k >= 2 else float("nan")
        ch  = ch_medoids(D, labels, medoids) if k >= 2 else float("nan")
        db  = db_medoids(D, labels, medoids) if k >= 2 else float("nan")
        results.append(dict(k=k, labels=labels, medoids=medoids, obj=obj, iters=iters,
                            sil=sil, ch=ch, db=db, backend=backend))
    results.sort(key=lambda r: (-np.nan_to_num(r["sil"], nan=-1e9),
                                -np.nan_to_num(r["ch"],  nan=-1e9),
                                 np.nan_to_num(r["db"],  nan=+1e9)))
    return results[0], results  # best, all


# ----------------------------- Run discovery -----------------------------

def find_runs(roots: List[str]) -> List[str]:
    runs = set()
    for root in roots:
        root = os.path.abspath(root)
        if not os.path.exists(root):
            continue
        for dirpath, dirnames, filenames in os.walk(root):
            if "wmfd" in dirnames:
                runs.add(dirpath)
    return sorted(runs)

def pick_matrix(run_dir: str) -> Optional[str]:
    cand_clean = os.path.join(run_dir, "wmfd", "trees_wmfd_matrix.csv")
    cand_noisy = os.path.join(run_dir, "wmfd", "trees_noisy_wmfd_matrix.csv")
    if os.path.exists(cand_clean):
        return cand_clean
    if os.path.exists(cand_noisy):
        return cand_noisy
    return None


# ----------------------------- Main -----------------------------

def main():
    cwd = os.path.abspath(os.getcwd())
    roots = [cwd, os.path.join(cwd, "A_noise_dense")]

    print("=== 1c auto K-Medoids (no-args) ===")
    print(f"Python : {'.'.join(map(str, __import__('sys').version_info[:3]))}")
    print("Roots:")
    for r in roots:
        print(f"  - {r}")

    runs = find_runs(roots)
    print(f"Runs found: {len(runs)}")
    if not runs:
        print("No runs found. Make sure your runs live under one of the roots above.")
        return

    summary_path = os.path.join(cwd, "kmedoids_summary.csv")
    with open(summary_path, "w", newline="", encoding="utf-8") as fsum:
        wsum = csv.writer(fsum)
        wsum.writerow([
            "run_dir","matrix_file","n","best_k","silhouette","calinski_harabasz","davies_bouldin",
            "objective","n_iter","backend","has_labels","ARI_bestK"
        ])

        for rdir in runs:
            mpath = pick_matrix(rdir)
            if not mpath:
                print(f"[skip] no matrix in {rdir}")
                continue

            try:
                ids, D = read_distance_matrix(mpath)
            except Exception as e:
                print(f"[error] read matrix failed in {rdir}: {e}")
                continue

            n = D.shape[0]
            K_list = list(range(2, max(3, min(12, n - 1)) + 1))  # 2..min(12, n-1)

            best, _all = choose_best_K(D, K_list, seed=0, max_iter=100)

            # Optional ARI if labels are present
            y_true = (read_true_labels(os.path.join(rdir, "trees.csv"), ids)
                      or read_true_labels(os.path.join(rdir, "labels.csv"), ids)
                      or read_true_labels(os.path.join(rdir, "wmfd", "labels.csv"), ids))
            ari_best = adjusted_rand_index(y_true, best["labels"]) if y_true is not None else None

            out_dir = os.path.join(rdir, "cluster_kmedoids", f"K={best['k']}")
            write_clusters(os.path.join(out_dir, "clusters_pred.csv"), ids, best["labels"])
            write_medoids(os.path.join(out_dir, "medoids.csv"), best["medoids"], ids)
            write_dist_to_medoids(os.path.join(out_dir, "dist_to_medoids.csv"), ids, best["medoids"], D)
            write_metrics(os.path.join(out_dir, "metrics.csv"), dict(
                k=best["k"], objective=f"{best['obj']:.6f}", n_iter=best["iters"], n_init=10, seed=0,
                ARI=(f"{ari_best:.6f}" if ari_best is not None else ""),
                silhouette=(f"{best['sil']:.6f}" if not math.isnan(best["sil"]) else ""),
                calinski_harabasz=(f"{best['ch']:.6f}" if not math.isnan(best["ch"]) else ""),
                davies_bouldin=(f"{best['db']:.6f}" if not math.isnan(best["db"]) else ""),
                backend=best["backend"]
            ))

            wsum.writerow([
                rdir, mpath, n, best["k"],
                (f"{best['sil']:.6f}" if not math.isnan(best["sil"]) else ""),
                (f"{best['ch']:.6f}" if not math.isnan(best["ch"]) else ""),
                (f"{best['db']:.6f}" if not math.isnan(best["db"]) else ""),
                f"{best['obj']:.6f}", best["iters"], best["backend"],
                ("yes" if y_true is not None else "no"),
                (f"{ari_best:.6f}" if ari_best is not None else "")
            ])

            sil_str = f"{best['sil']:.4f}" if not math.isnan(best["sil"]) else "nan"
            print(f"[OK] {os.path.basename(rdir)}  n={n}  bestK={best['k']}  sil={sil_str}")

    print(f"\nSummary written to: {summary_path}")


if __name__ == "__main__":
    main()