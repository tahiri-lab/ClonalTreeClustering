# -- coding: utf-8 --
"""
wmfd_from_any.py — ZERO ARGUMENTS
Automatically reads an input file,
extracts only Newick strings (no other columns required),
computes the WMFD distance matrix, and writes:
  - trees_wmfd_matrix.csv
  - trees_wmfd_pairs.csv

Dependency: ete3  (pip install ete3)
"""

import os
import re
import csv
from typing import List, Tuple, Dict, Set
from ete3 import Tree

# ----------- Newick detection helpers -----------

def _looks_like_newick(s: str) -> bool:
    s = s.strip().strip('"').strip("'")
    return "(" in s and ")" in s and s.rstrip().endswith(";")

def _extract_first_newick(s: str) -> str | None:
    s0 = s.strip().strip('"').strip("'")
    if not s0:
        return None
    m = re.search(r"\(.*?\);", s0, flags=re.DOTALL)
    return m.group(0).replace("\n", "").replace("\r", "") if m else None

def read_newicks_any(path: str) -> List[str]:
    """Try line-by-line, then CSV cells, then a global regex extraction."""
    newicks: List[str] = []

    # 1) line-by-line
    try:
        with open(path, "r", encoding="utf-8", newline="") as f:
            for raw in f:
                if not raw.strip():
                    continue
                if _looks_like_newick(raw):
                    cand = _extract_first_newick(raw)
                    if not cand:
                        continue
                    try:
                        _ = Tree(cand, format=1)
                        newicks.append(cand)
                    except Exception:
                        pass
    except Exception:
        pass
    if newicks:
        return newicks

    # 2) parse CSV cells
    try:
        with open(path, "r", encoding="utf-8", newline="") as f:
            rdr = csv.reader(f)
            for row in rdr:
                for cell in row:
                    if _looks_like_newick(cell):
                        cand = _extract_first_newick(cell)
                        if not cand:
                            continue
                        try:
                            _ = Tree(cand, format=1)
                            newicks.append(cand)
                        except Exception:
                            pass
    except Exception:
        pass
    if newicks:
        return newicks

    # 3) global extraction
    try:
        with open(path, "r", encoding="utf-8", newline="") as f:
            txt = f.read()
        for cand in re.findall(r"\(.*?\);", txt, flags=re.DOTALL):
            c = cand.replace("\n", "").replace("\r", "").strip()
            if not c:
                continue
            try:
                _ = Tree(c, format=1)
                newicks.append(c)
            except Exception:
                pass
    except Exception:
        pass

    return newicks

# ----------------------------- WMFD -----------------------------

def _parse_leaf_value(name: str) -> float:
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

def wmfd_pair(A, B, L1=0.30, L2=0.20, L3=0.25, L4=0.15, L5=0.10) -> float:
    (_t1, bl1, h1, w1, d1, Ls1, S1) = A
    (_t2, bl2, h2, w2, d2, Ls2, S2) = B

    tot = L1 + L2 + L3 + L4 + L5
    if tot <= 0:
        L1, L2, L3, L4, L5 = 0.30, 0.20, 0.25, 0.15, 0.10
        tot = 1.0
    L1, L2, L3, L4, L5 = (L1/tot, L2/tot, L3/tot, L4/tot, L5/tot)

    unionL = Ls1 | Ls2
    interL = Ls1 & Ls2
    TN = len(unionL); CN = len(interL)
    P = 1.0 - (CN / TN) if TN > 0 else 0.0

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

    CN = len(Ls1 & Ls2)
    TN = len(Ls1 | Ls2)
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

    S1 = A[6]; S2 = B[6]
    unionS = S1 | S2
    interS = S1 & S2
    HD = 1.0 - (len(interS)/len(unionS)) if len(unionS) > 0 else 0.0

    P = 1.0 - (CN / TN) if TN > 0 else 0.0
    return float(P * WND_uncommon + WND_common + L5 * HD)

# ----------------------------- Writers -----------------------------

def write_outputs(stem: str, D):
    n = len(D)
    out_matrix = f"{stem}_wmfd_matrix.csv"
    out_pairs  = f"{stem}_wmfd_pairs.csv"

    # --- Matrix file ---
    with open(out_matrix, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        for i in range(n):
            w.writerow([f"{D[i][j]:.8f}" for j in range(n)])

    # --- Pairs file ---
    with open(out_pairs, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        for i in range(n):
            for j in range(i+1, n):
                w.writerow([f"{D[i][j]:.8f}"])

    return out_matrix, out_pairs

# ----------------------------- Main -----------------------------

if __name__ == "__main__":
    # Hard-coded input candidates
    candidates = ["trees_4_16_17_8_30.txt"]
    infile = None
    for c in candidates:
        p = os.path.abspath(c)
        if os.path.exists(p):
            infile = p
            break

    if not infile:
        print("No input file found (looked for: trees_4_16_17_8_30.txt).")
        raise SystemExit(1)

    print(f"[LOAD] {infile}")
    newicks = read_newicks_any(infile)
    if not newicks:
        print("[ERR] No valid Newick strings detected in the file.")
        raise SystemExit(1)

    print(f"[OK] {len(newicks)} Newick trees detected.")
    feats = [precompute_features(nw) for nw in newicks]

    n = len(feats)
    D = [[0.0]*n for _ in range(n)]

    print("[WMFD] computing pairwise distances…")
    for i in range(n):
        if i % 5 == 0:
            print(f"  {i}/{n}")
        for j in range(i+1, n):
            d = wmfd_pair(feats[i], feats[j])
            D[i][j] = D[j][i] = float(d)
    print(f"  {n}/{n}")

    stem = os.path.splitext(infile)[0]  # e.g., "trees"
    out_matrix, out_pairs = write_outputs(stem, D)
    print(f"[DONE] Matrix : {out_matrix}")
    print(f"[DONE] Pairs  : {out_pairs}")