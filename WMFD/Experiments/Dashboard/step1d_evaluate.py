# -*- coding: utf-8 -*-
"""
Step 1d) Evaluate clustering: ARI, contingency, purity, mapping (majority), errors list.

Entries:
- clusters_pred.csv: id,pred_cluster
- labels.csv: id, true_cluster

Outputs:
- ari_report.txt: ARI + purity + summary
- contingency.csv: matrix (pred x true)
- mapping_majority.csv: pred_cluster -> true_cluster_majority
- misassigned.csv: id, true, pred (only if error)

Use:
python step1d_evaluate.py --pred clusters_pred.csv --true labels.csv
"""

import csv
import argparse
from collections import defaultdict

def _comb2(x): return x*(x-1)/2.0

def adjusted_rand_index(y_true, y_pred):
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
    n_pairs = _comb2(len(y_true))
    expected = (ai_sum * bj_sum) / n_pairs if n_pairs else 0.0
    max_index = 0.5*(ai_sum + bj_sum)
    denom = max_index - expected
    if denom == 0:
        return 1.0 if nij_sum == max_index else 0.0
    return (nij_sum - expected) / denom

def read_pred(path):
    ids, ypred = [], []
    with open(path, newline="", encoding="utf-8") as f:
        r = csv.DictReader(f)
        for row in r:
            ids.append(row["id"])
            ypred.append(int(row["pred_cluster"]))
    return ids, ypred

def read_true(path, want_ids):
    m = {}
    with open(path, newline="", encoding="utf-8") as f:
        r = csv.DictReader(f)
        for row in r:
            m[row["id"]] = int(row["true_cluster"])
    return [m[i] for i in want_ids]

def contingency(ids, y_true, y_pred):
    # cluster keys
    C_true = sorted(set(y_true))
    C_pred = sorted(set(y_pred))
    idx_true = {c:i for i,c in enumerate(C_true)}
    idx_pred = {c:i for i,c in enumerate(C_pred)}
    M = [[0]*len(C_true) for _ in range(len(C_pred))]
    for t,p in zip(y_true, y_pred):
        M[idx_pred[p]][idx_true[t]] += 1
    return C_pred, C_true, M

def majority_mapping(C_pred, C_true, M):
    mapping = {}
    for ip, p in enumerate(C_pred):
        # max column
        j = max(range(len(C_true)), key=lambda j: M[ip][j])
        mapping[p] = C_true[j]
    return mapping  #dict pred -> true_majority

def purity(M):
    tot = sum(sum(row) for row in M)
    maj_sum = sum(max(row) for row in M)
    return maj_sum / tot if tot else 0.0

def write_contingency(path, C_pred, C_true, M):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["pred\\true"] + [str(c) for c in C_true])
        for ip, p in enumerate(C_pred):
            w.writerow([p] + M[ip])

def write_mapping(path, mapping):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["pred_cluster","true_majority"])
        for p in sorted(mapping.keys()):
            w.writerow([p, mapping[p]])

def write_errors(path, ids, y_true, y_pred):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["id","true","pred"])
        for i,(t,p) in enumerate(zip(y_true, y_pred)):
            if t != p:
                w.writerow([ids[i], t, p])

def main():
    ap = argparse.ArgumentParser(description="Evaluate clustering (ARI, contingency, purity).")
    ap.add_argument("--pred", required=True, help="clusters_pred.csv")
    ap.add_argument("--true", required=True, help="labels.csv")
    ap.add_argument("--out_report", default="ari_report.txt")
    ap.add_argument("--out_cont", default="contingency.csv")
    ap.add_argument("--out_map",  default="mapping_majority.csv")
    ap.add_argument("--out_err",  default="misassigned.csv")
    args = ap.parse_args()

    ids, y_pred = read_pred(args.pred)
    y_true = read_true(args.true, ids)

    ari = adjusted_rand_index(y_true, y_pred)
    C_pred, C_true, M = contingency(ids, y_true, y_pred)
    pur = purity(M)
    mapping = majority_mapping(C_pred, C_true, M)

    write_contingency(args.out_cont, C_pred, C_true, M)
    write_mapping(args.out_map, mapping)
    write_errors(args.out_err, ids, y_true, y_pred)

    with open(args.out_report, "w", encoding="utf-8") as f:
        f.write(f"ARI: {ari:.6f}\n")
        f.write(f"Purity (majority): {pur:.6f}\n")
        f.write(f"Pred clusters: {C_pred}\n")
        f.write(f"True clusters: {C_true}\n")
        f.write("Majority mapping (pred -> true):\n")
        for p in sorted(mapping.keys()):
            f.write(f"  {p} -> {mapping[p]}\n")
    print(f"Done. ARI={ari:.6f}  Purity={pur:.6f}")

if __name__ == "__main__":
    main()
