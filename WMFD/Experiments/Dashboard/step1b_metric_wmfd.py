# -- coding: utf-8 --
"""
STEP 1b — WMFD metric on all runs (auto-discovery, recursive)

What it does
------------
- Recursively scans run folders starting from:
    1) current working directory,
    2) "<cwd>/A_noise_dense" if it exists.
- For each run folder, if it contains trees.csv and/or trees_noisy.csv,
  it computes the WMFD pairwise distance matrix.
- Writes outputs next to each input, under a wmfd/ subfolder:
    - <name>_wmfd_matrix.csv  (NxN matrix with headers)
    - <name>_wmfd_pairs.csv   (upper-triangle i<j)
    - <name>_labels.csv       (id,true_cluster)

Inputs (per file)
-----------------
CSV with columns (case/space/BOM tolerant):
  - cluster | true_cluster | c
  - tree_id | id | treeid | t
  - newick  | nwk | tree | newick_str

Outputs (per input file)
------------------------
{run_dir}/wmfd/
  - trees_wmfd_matrix.csv
  - trees_wmfd_pairs.csv
  - trees_labels.csv
and/or
  - trees_noisy_wmfd_matrix.csv
  - trees_noisy_wmfd_pairs.csv
  - trees_noisy_labels.csv

Notes
-----
- WMFD uses pairwise min-max normalization per channel (BL/H/W/D).
- If you re-run, files are overwritten by default.
- No PowerShell args needed; just python step1b_metric_wmfd.py.
"""

import os
import csv
import sys
from typing import Dict, Set, List, Tuple
from ete3 import Tree

# ----------------------------- IO utils -----------------------------

def _norm_keys(fieldnames):
    """Original header key -> normalized (strip, lower, BOM removed)."""
    out = {}
    for k in fieldnames or []:
        out[k] = (k or "").strip().lower().replace("\ufeff", "")
    return out

def _pick_key(keys_map, wanted):
    """Return first original key whose normalized version is in 'wanted'."""
    wanted = set(wanted)
    for orig, norm in keys_map.items():
        if norm in wanted:
            return orig
    return None

def read_trees_csv(in_csv: str):
    """Read a trees CSV and validate required columns + Newick parsing."""
    if not os.path.exists(in_csv):
        raise FileNotFoundError(f"Missing file: {in_csv}")

    ids: List[str] = []
    clusters: List[int] = []
    newicks: List[str] = []

    with open(in_csv, newline="", encoding="utf-8") as f:
        r = csv.DictReader(f)
        if not r.fieldnames:
            raise ValueError(f"Missing header in: {in_csv}")

        keys = _norm_keys(r.fieldnames)
        c_key  = _pick_key(keys, {"cluster", "true_cluster", "c"})
        t_key  = _pick_key(keys, {"tree_id", "id", "treeid", "t"})
        nw_key = _pick_key(keys, {"newick", "nwk", "tree", "newick_str"})

        if not (c_key and t_key and nw_key):
            raise KeyError(
                f"Required headers not found. Present={r.fieldnames}. "
                "Accepted: cluster|true_cluster, tree_id|id, newick|nwk|tree."
            )

        for row in r:
            cs = str(row[c_key]).strip()
            ts = str(row[t_key]).strip()
            nw = str(row[nw_key]).strip()
            if not cs or not ts or not nw:
                continue
            try:
                c = int(float(cs))
                t = int(float(ts))
            except ValueError:
                raise ValueError(f"Invalid row (cluster='{cs}', tree_id='{ts}') in {in_csv}")

            # Validate Newick quickly
            try:
                _ = Tree(nw, format=1)
            except Exception as e:
                raise ValueError(f"Newick parse error for id={c}.{t} in {in_csv}: {e}")

            ids.append(f"{c}.{t}")
            clusters.append(c)
            newicks.append(nw)

    if not ids:
        raise ValueError(f"No readable trees in {in_csv}")
    return ids, clusters, newicks

# ----------------------------- WMFD helpers -----------------------------

def _parse_leaf_value(name: str) -> float:
    # "L12@3.7" -> 3.7 ; "L12" -> 1.0 (default)
    try:
        if "@" in name:
            return float(name.split("@", 1)[1])
    except Exception:
        pass
    return 1.0

def _parent_degree(leaf) -> int:
    p = leaf.up
    return len(p.children) if p is not None else 1

def precompute_features(nwk: str):
    """
    Return (tree, leaf_BL, leaf_H, leaf_W, leaf_D, leaf_set, splits_set).
    BL/H/W/D are raw here; pairwise min-max normalization is done in wmfd_pair.
    """
    t = Tree(nwk, format=1)
    root = t.get_tree_root()

    leaf_BL: Dict[str, float] = {}
    leaf_H: Dict[str, float]  = {}
    leaf_W: Dict[str, float]  = {}
    leaf_D: Dict[str, float]  = {}

    for leaf in t.iter_leaves():
        name = str(leaf.name)
        leaf_BL[name] = float(leaf.dist or 0.0)
        leaf_H[name]  = float(root.get_distance(leaf) or 0.0)
        leaf_W[name]  = float(_parse_leaf_value(name))
        leaf_D[name]  = float(_parent_degree(leaf))

    leaf_set: Set[str] = set(leaf_BL.keys())

    splits: Set[frozenset] = set()
    for node in t.traverse("postorder"):
        if node.is_leaf():
            continue
        clade = frozenset(l.name for l in node.iter_leaves())
        if 1 < len(clade) < len(leaf_set):
            splits.add(clade)

    return t, leaf_BL, leaf_H, leaf_W, leaf_D, leaf_set, splits

def _minmax_pair_norm(v1: float, v2: float, mn: float, mx: float):
    if mx > mn:
        return ((v1 - mn) / (mx - mn), (v2 - mn) / (mx - mn))
    return (0.0, 0.0)

def wmfd_pair(A, B, L1: float, L2: float, L3: float, L4: float, L5: float) -> float:
    """
    WMFD = P * WND_uncommon + WND_common + L5 * HD
    with pairwise min-max normalization per channel (BL/H/W/D).
    """
    (_t1, bl1, h1, w1, d1, Ls1, S1) = A
    (_t2, bl2, h2, w2, d2, Ls2, S2) = B

    # Normalize weights to sum to 1
    tot = L1 + L2 + L3 + L4 + L5
    if tot <= 0:
        L1, L2, L3, L4, L5 = 0.30, 0.20, 0.25, 0.15, 0.10
        tot = 1.0
    L1, L2, L3, L4, L5 = (L1/tot, L2/tot, L3/tot, L4/tot, L5/tot)

    unionL = Ls1 | Ls2
    interL = Ls1 & Ls2
    TN = len(unionL); CN = len(interL)
    P = 1.0 - (CN / TN) if TN > 0 else 0.0

    # min-max bounds over the union
    bl_min = min([bl1.get(x,0.0) for x in unionL] + [bl2.get(x,0.0) for x in unionL], default=0.0)
    bl_max = max([bl1.get(x,0.0) for x in unionL] + [bl2.get(x,0.0) for x in unionL], default=0.0)
    h_min  = min([h1.get(x,0.0)  for x in unionL] + [h2.get(x,0.0)  for x in unionL], default=0.0)
    h_max  = max([h1.get(x,0.0)  for x in unionL] + [h2.get(x,0.0)  for x in unionL], default=0.0)
    w_min  = min([w1.get(x,0.0)  for x in unionL] + [w2.get(x,0.0)  for x in unionL], default=0.0)
    w_max  = max([w1.get(x,0.0)  for x in unionL] + [w2.get(x,0.0)  for x in unionL], default=0.0)
    d_min  = min([d1.get(x,0.0)  for x in unionL] + [d2.get(x,0.0)  for x in unionL], default=0.0)
    d_max  = max([d1.get(x,0.0)  for x in unionL] + [d2.get(x,0.0)  for x in unionL], default=0.0)

    com_BL = com_H = com_W = com_D = 0.0
    unc_BL = unc_H = unc_W = unc_D = 0.0

    for leaf in unionL:
        bl_a, bl_b = _minmax_pair_norm(bl1.get(leaf,0.0), bl2.get(leaf,0.0), bl_min, bl_max)
        h_a,  h_b  = _minmax_pair_norm(h1.get(leaf,0.0),  h2.get(leaf,0.0),  h_min,  h_max)
        w_a,  w_b  = _minmax_pair_norm(w1.get(leaf,0.0),  w2.get(leaf,0.0),  w_min,  w_max)
        d_a,  d_b  = _minmax_pair_norm(d1.get(leaf,0.0),  d2.get(leaf,0.0),  d_min,  d_max)

        diff_BL = abs(bl_a - bl_b)
        diff_H  = abs(h_a  - h_b)
        diff_W  = abs(w_a  - w_b)
        diff_D  = abs(d_a  - d_b)

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

    # HD on splits (Jaccard)
    unionS = S1 | S2
    interS = S1 & S2
    HD = 1.0 - (len(interS)/len(unionS)) if len(unionS) > 0 else 0.0

    return float(P * WND_uncommon + WND_common + L5 * HD)

# ----------------------------- Writers -----------------------------

def write_outputs(ids, clusters, D, out_matrix, out_pairs, out_labels):
    n = len(ids)
    os.makedirs(os.path.dirname(out_matrix), exist_ok=True)

    with open(out_matrix, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["id"] + ids)
        for i in range(n):
            w.writerow([ids[i]] + [("{:.8f}".format(D[i][j])) for j in range(n)])

    with open(out_pairs, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["id_i","id_j","dist"])
        for i in range(n):
            for j in range(i+1, n):
                w.writerow([ids[i], ids[j], "{:.8f}".format(D[i][j])])

    with open(out_labels, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["id","true_cluster"])
        for i, lab in enumerate(clusters):
            w.writerow([ids[i], lab])

# ----------------------------- Core runner --------------------------

def process_one_csv(in_csv: str, L1=0.30, L2=0.20, L3=0.25, L4=0.15, L5=0.10):
    ids, clusters, newicks = read_trees_csv(in_csv)
    feats = [precompute_features(nwk) for nwk in newicks]

    n = len(ids)
    D = [[0.0]*n for _ in range(n)]
    base = os.path.basename(in_csv)
    print(f"[WMFD] {base}  0/{n}")
    for i in range(n):
        if (i % 5) == 0 and i > 0:
            print(f"[WMFD] {base}  {i}/{n}")
        for j in range(i+1, n):
            d = wmfd_pair(feats[i], feats[j], L1, L2, L3, L4, L5)
            D[i][j] = D[j][i] = float(d)
    print(f"[WMFD] {base}  {n}/{n}")

    run_dir   = os.path.dirname(in_csv)
    out_dir   = os.path.join(run_dir, "wmfd")
    stem      = os.path.splitext(os.path.basename(in_csv))[0]  # "trees" or "trees_noisy"
    out_matrix = os.path.join(out_dir, f"{stem}_wmfd_matrix.csv")
    out_pairs  = os.path.join(out_dir, f"{stem}_wmfd_pairs.csv")
    out_labels = os.path.join(out_dir, f"{stem}_labels.csv")

    write_outputs(ids, clusters, D, out_matrix, out_pairs, out_labels)
    print(f"[OK] WMFD → {out_matrix}")

# ----------------------------- Discovery ----------------------------

def find_roots() -> List[str]:
    roots = []
    cwd = os.getcwd()
    roots.append(cwd)
    dense = os.path.join(cwd, "A_noise_dense")
    if os.path.isdir(dense):
        roots.append(dense)
    # also allow walking up one level (useful if script sits in a subdir)
    up1 = os.path.dirname(cwd)
    if up1 and os.path.isdir(up1):
        roots.append(up1)
        up1_dense = os.path.join(up1, "A_noise_dense")
        if os.path.isdir(up1_dense):
            roots.append(up1_dense)
    # deduplicate while preserving order
    seen = set()
    uniq = []
    for r in roots:
        if r not in seen:
            uniq.append(r); seen.add(r)
    return uniq

def discover_input_csvs(roots: List[str]) -> List[str]:
    """
    Recursively walk roots and collect paths to 'trees.csv' and 'trees_noisy.csv'.
    """
    targets = []
    for root in roots:
        for dirpath, _dirnames, filenames in os.walk(root):
            if "wmfd" in dirpath.replace("\\", "/").split("/"):
                # skip wmfd output folders
                continue
            files = set(filenames)
            for name in ("trees.csv", "trees_noisy.csv"):
                if name in files:
                    targets.append(os.path.join(dirpath, name))
    # remove duplicates, keep order
    seen = set(); out = []
    for p in targets:
        if p not in seen:
            out.append(p); seen.add(p)
    return out

# ----------------------------- Main --------------------------------

def main():
    print("=== STEP 1b (WMFD) — recursive auto-discovery ===")
    print("Python :", sys.version)
    print("CWD    :", os.getcwd())

    roots = find_roots()
    print("Search roots:")
    for r in roots:
        print("  -", r)

    csvs = discover_input_csvs(roots)
    if not csvs:
        print("\nNo CSV found. Put your runs under one of the roots above "
              "(ideally ...\\A_noise_dense) and re-run.")
        return

    # Process "trees.csv" before "trees_noisy.csv" for each run directory
    # by ordering the list: normal first, noisy second when same folder.
    def sort_key(p):
        base = os.path.basename(p)
        return (os.path.dirname(p), 0 if base == "trees.csv" else 1, base)
    csvs.sort(key=sort_key)

    print("\nInputs to process:")
    for p in csvs:
        print("  ·", p)

    for path in csvs:
        try:
            process_one_csv(path)
        except Exception as e:
            print(f"[WARN] {path}: {e}")

if __name__ == "__main__":
    main()