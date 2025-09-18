# step1d_evaluate.py
# -- coding: utf-8 --
"""
Step 1d) Evaluate clustering: ARI, contingency, purity, majority mapping, error list.

Two modes:
1) Single run:
   python step1d_evaluate.py --pred <.../clusters_pred.csv> --true <.../labels.csv>

2) Batch (recommended): scan all runs under root (default: ./A_noise_dense)
   python step1d_evaluate.py --root ./A_noise_dense

Inputs expected per run
----------------------
- clusters_pred.csv: columns {id, pred_cluster}
  (either directly in the run directory, or in cluster_kmedoids/K=*/clusters_pred.csv)
- labels.csv or trees.csv: columns {id, true_cluster} (tolerant to header variants)

Per-run outputs (written under <run>/eval/)
-------------------------------------------
- ari_report.txt              : text summary with ARI, purity, mapping
- contingency.csv             : contingency table (pred x true)
- mapping_majority.csv        : mapping pred_cluster -> true_majority
- misassigned.csv             : only misclassified items

Global summary
--------------
- evaluation_summary.csv      : one line per run (path, ARI, purity, #items, K if known)

Notes
-----
- Header parsing is robust (case/whitespace/BOM tolerant).
- If multiple predicted clusterings exist (e.g. several K under cluster_kmedoids), we try:
    (a) pick the one with the best silhouette in its metrics.csv if available,
    else (b) pick the most recently modified clusters_pred.csv.
"""

import os
import csv
import glob
import math
import argparse
from collections import defaultdict
from typing import Dict, List, Tuple, Optional


# ----------------------------- small utils -----------------------------

def _comb2(x: int) -> float:
    return x * (x - 1) / 2.0

def _norm_keys(fieldnames):
    out = {}
    for k in fieldnames or []:
        out[k] = (k or "").strip().lower().replace("\ufeff", "")
    return out

def _to_int(x: str) -> int:
    s = (x or "").strip().replace(",", ".")
    return int(float(s))

def _ensure_dir(p: str):
    os.makedirs(p, exist_ok=True)


# ----------------------------- ARI & helpers -----------------------------

def adjusted_rand_index(y_true: List[int], y_pred: List[int]) -> float:
    n = len(y_true)
    if n == 0:
        return float("nan")
    by_t = defaultdict(list)
    by_p = defaultdict(list)
    for i, (t, p) in enumerate(zip(y_true, y_pred)):
        by_t[t].append(i)
        by_p[p].append(i)
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
    max_index = 0.5 * (ai_sum + bj_sum)
    denom = max_index - expected
    if denom == 0:
        return 1.0 if nij_sum == max_index else 0.0
    return (nij_sum - expected) / denom


def contingency(ids: List[str], y_true: List[int], y_pred: List[int]):
    C_true = sorted(set(y_true))
    C_pred = sorted(set(y_pred))
    idx_true = {c: i for i, c in enumerate(C_true)}
    idx_pred = {c: i for i, c in enumerate(C_pred)}
    M = [[0] * len(C_true) for _ in range(len(C_pred))]
    for t, p in zip(y_true, y_pred):
        M[idx_pred[p]][idx_true[t]] += 1
    return C_pred, C_true, M

def purity(M: List[List[int]]) -> float:
    tot = sum(sum(row) for row in M)
    if tot == 0:
        return float("nan")
    maj = sum(max(row) for row in M) if M else 0
    return maj / tot

def majority_mapping(C_pred: List[int], C_true: List[int], M: List[List[int]]) -> Dict[int, int]:
    mapping = {}
    for ip, p in enumerate(C_pred):
        j = max(range(len(C_true)), key=lambda j: M[ip][j]) if C_true else None
        mapping[p] = C_true[j] if j is not None else None
    return mapping


# ----------------------------- I/O: read predictors & labels -----------------------------

def read_pred_clusters(pred_csv: str) -> Tuple[List[str], List[int]]:
    with open(pred_csv, newline="", encoding="utf-8") as f:
        rdr = csv.DictReader(f)
        if not rdr.fieldnames:
            raise ValueError(f"Missing header in {pred_csv}")
        keys = _norm_keys(rdr.fieldnames)
        id_key = next((k for k, v in keys.items() if v == "id"), None)
        pred_key = next((k for k, v in keys.items()
                         if v in ("pred_cluster", "pred", "cluster", "label", "y")), None)
        if not id_key or not pred_key:
            raise KeyError(f"{pred_csv}: need columns 'id' and 'pred_cluster' (tolerant). Found {rdr.fieldnames}")
        ids, ypred = [], []
        for row in rdr:
            ids.append(str(row[id_key]).strip())
            ypred.append(_to_int(row[pred_key]))
    return ids, ypred

def read_true_labels(true_csv: str, ids: List[str]) -> List[int]:
    with open(true_csv, newline="", encoding="utf-8") as f:
        rdr = csv.DictReader(f)
        if not rdr.fieldnames:
            raise ValueError(f"Missing header in {true_csv}")
        keys = _norm_keys(rdr.fieldnames)
        id_key = next((k for k, v in keys.items() if v == "id"), None)
        true_key = next((k for k, v in keys.items()
                         if v in ("true_cluster", "true", "label", "cluster", "y")), None)
        if not id_key or not true_key:
            raise KeyError(f"{true_csv}: need columns 'id' and 'true_cluster' (tolerant). Found {rdr.fieldnames}")
        mapping = {}
        for row in rdr:
            mapping[str(row[id_key]).strip()] = _to_int(row[true_key])
    # order to match ids from pred file
    try:
        return [mapping[i] for i in ids]
    except KeyError as e:
        missing = str(e).strip("'")
        raise KeyError(f"{true_csv}: missing label for id '{missing}'")

# Try to find a “true labels” file for a run directory.
def auto_pick_true_csv(run_dir: str) -> Optional[str]:
    # prefer explicit labels.csv at run root
    cand = os.path.join(run_dir, "labels.csv")
    if os.path.exists(cand):
        return cand
    # sometimes under wmfd/
    cand = os.path.join(run_dir, "wmfd", "labels.csv")
    if os.path.exists(cand):
        return cand
    # fallback: trees.csv with column 'cluster'
    cand = os.path.join(run_dir, "trees.csv")
    if os.path.exists(cand):
        return cand
    # even noisier variant
    cand = os.path.join(run_dir, "trees_noisy.csv")
    if os.path.exists(cand):
        return cand
    return None

# Find the best clusters_pred.csv for a run.
def auto_pick_pred_csv(run_dir: str) -> Optional[str]:
    # 1) direct file at run root
    direct = os.path.join(run_dir, "clusters_pred.csv")
    if os.path.exists(direct):
        return direct
    # 2) search under cluster_kmedoids/K=*/clusters_pred.csv
    pattern = os.path.join(run_dir, "cluster_kmedoids", "K=*",
                           "clusters_pred.csv")
    cands = glob.glob(pattern)
    if not cands:
        return None
    # If metrics.csv exists alongside, pick best by silhouette; else pick most recent.
    scored = []
    for p in cands:
        mpath = os.path.join(os.path.dirname(p), "metrics.csv")
        sil = None
        if os.path.exists(mpath):
            try:
                with open(mpath, newline="", encoding="utf-8") as f:
                    rows = list(csv.DictReader(f))
                    if rows:
                        s = (rows[0].get("silhouette") or "").strip()
                        sil = float(s) if s else None
            except Exception:
                sil = None
        scored.append((p, sil, os.path.getmtime(p)))
    # Sort: higher silhouette first (None last), then latest mtime
    scored.sort(key=lambda t: ((-t[1]) if t[1] is not None else math.inf, -t[2]))
    return scored[0][0]


# ----------------------------- writers -----------------------------

def write_contingency(path: str, C_pred: List[int], C_true: List[int], M: List[List[int]]) -> None:
    _ensure_dir(os.path.dirname(path))
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["pred\\true"] + [str(c) for c in C_true])
        for ip, p in enumerate(C_pred):
            w.writerow([p] + M[ip])

def write_mapping(path: str, mapping: Dict[int, int]) -> None:
    _ensure_dir(os.path.dirname(path))
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["pred_cluster", "true_majority"])
        for p in sorted(mapping.keys()):
            w.writerow([p, mapping[p]])

def write_errors(path: str, ids: List[str], y_true: List[int], y_pred: List[int]) -> None:
    _ensure_dir(os.path.dirname(path))
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["id", "true", "pred"])
        for i, (t, p) in enumerate(zip(y_true, y_pred)):
            if t != p:
                w.writerow([ids[i], t, p])

def write_report(path: str, ari: float, pur: float, mapping: Dict[int, int],
                 C_pred: List[int], C_true: List[int]) -> None:
    _ensure_dir(os.path.dirname(path))
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"ARI: {ari:.6f}\n")
        f.write(f"Purity (majority): {pur:.6f}\n")
        f.write(f"Pred clusters: {C_pred}\n")
        f.write(f"True clusters: {C_true}\n")
        f.write("Majority mapping (pred -> true):\n")
        for p in sorted(mapping.keys()):
            f.write(f"  {p} -> {mapping[p]}\n")


# ----------------------------- single-run evaluation -----------------------------

def evaluate_one(pred_csv: str, true_csv: str, out_dir: str) -> Tuple[int, float, float]:
    ids, y_pred = read_pred_clusters(pred_csv)
    y_true = read_true_labels(true_csv, ids)

    C_pred, C_true, M = contingency(ids, y_true, y_pred)
    ari = adjusted_rand_index(y_true, y_pred)
    pur = purity(M)
    mapping = majority_mapping(C_pred, C_true, M)

    # write outputs
    write_contingency(os.path.join(out_dir, "contingency.csv"), C_pred, C_true, M)
    write_mapping(os.path.join(out_dir, "mapping_majority.csv"), mapping)
    write_errors(os.path.join(out_dir, "misassigned.csv"), ids, y_true, y_pred)
    write_report(os.path.join(out_dir, "ari_report.txt"), ari, pur, mapping, C_pred, C_true)

    return len(ids), ari, pur


# ----------------------------- batch driver -----------------------------

def scan_runs(root: str) -> List[str]:
    runs = []
    for dirpath, dirnames, filenames in os.walk(root):
        # A run directory should contain either clusters_pred.csv or a cluster_kmedoids subfolder
        if ("clusters_pred.csv" in filenames) or ("cluster_kmedoids" in dirnames):
            runs.append(dirpath)
    return sorted(set(runs))

def main():
    ap = argparse.ArgumentParser(description="1d) Evaluate clustering (ARI, purity, contingency, mapping).")
    ap.add_argument("--pred", type=str, help="clusters_pred.csv (single run mode)")
    ap.add_argument("--true", type=str, help="labels.csv or trees.csv (single run mode)")
    ap.add_argument("--root", type=str, default=None, help="Scan all runs under this folder (batch mode)")
    ap.add_argument("--summary", type=str, default="evaluation_summary.csv",
                    help="Global CSV summary when using --root")
    args = ap.parse_args()

    # Mode 1: explicit files
    if args.pred and args.true:
        out_dir = os.path.join(os.path.dirname(os.path.abspath(args.pred)), "eval")
        n, ari, pur = evaluate_one(args.pred, args.true, out_dir)
        print(f"Done. n={n}  ARI={ari:.6f}  Purity={pur:.6f}  → {out_dir}")
        return

    # Mode 2: batch
    root = args.root or os.path.join(".", "A_noise_dense")
    root = os.path.abspath(root)
    if not os.path.exists(root):
        print(f"Root folder not found: {root}")
        return

    runs = scan_runs(root)
    print(f"=== 1d) Batch evaluation ===")
    print(f"Root: {root}")
    print(f"Runs found: {len(runs)}")
    if not runs:
        return

    summary_path = os.path.join(os.getcwd(), args.summary)
    with open(summary_path, "w", newline="", encoding="utf-8") as fsum:
        w = csv.writer(fsum)
        w.writerow(["run_dir", "pred_file", "true_file", "n", "ARI", "Purity"])

        for rdir in runs:
            pred_csv = auto_pick_pred_csv(rdir)
            true_csv = auto_pick_true_csv(rdir)
            if not pred_csv or not true_csv:
                print(f"[skip] {os.path.basename(rdir)}: missing files "
                      f"(pred={bool(pred_csv)} true={bool(true_csv)})")
                continue
            try:
                out_dir = os.path.join(rdir, "eval")
                n, ari, pur = evaluate_one(pred_csv, true_csv, out_dir)
                w.writerow([rdir, pred_csv, true_csv, n, f"{ari:.6f}", f"{pur:.6f}"])
                print(f"[OK] {os.path.basename(rdir)}  n={n}  ARI={ari:.4f}  Purity={pur:.4f}")
            except Exception as e:
                print(f"[error] {os.path.basename(rdir)}: {e}")

    print(f"\nSummary written to: {summary_path}")


if __name__ == "__main__":
    main()