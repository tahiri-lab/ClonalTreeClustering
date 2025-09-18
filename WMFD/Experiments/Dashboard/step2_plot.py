# step2_plot.py
# -- coding: utf-8 --
"""
Dashboard (a)…(f) à partir de results_for_plot.csv produit par make_results_for_plot.py

Entrée (CSV) — colonnes attendues :
  metric,k,L,n,noise,noise_mode,plevel,repeat,ARI,objective,coherence_measured,seed

Panneaux :
 (a) ARI vs noise par K
 (b) ARI vs L à noise_ref
 (c) ARI vs n à noise_ref
 (d) ARI vs plevel (ou cohérence mesurée si dispo + --use_measured_coherence)
 (e) ARI vs K à noise_ref
 (f) ΔARI (noise_low − noise_high) par K
"""

import argparse, csv, math
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt

def safe_float(x):
    if x is None: return float("nan")
    xs = str(x).strip().replace(",", ".")
    if xs == "" or xs.lower() == "nan": return float("nan")
    return float(xs)

def safe_int(x):
    try: return int(float(x))
    except: return 0

def load_rows(path, metric_filter=None, noise_mode_filter=None, plevel_filter=None):
    rows = []
    with open(path, newline="", encoding="utf-8") as f:
        rdr = csv.DictReader(f)
        for r in rdr:
            metric = (r.get("metric") or "").strip().lower()
            k = safe_int(r.get("k"))
            L = safe_int(r.get("L"))
            n = safe_int(r.get("n"))
            noise = safe_float(r.get("noise"))
            plevel = safe_float(r.get("plevel"))
            noise_mode = (r.get("noise_mode") or "").strip().lower()
            repeat = safe_int(r.get("repeat"))
            ARI = safe_float(r.get("ARI"))
            coh = safe_float(r.get("coherence_measured"))
            if metric_filter and metric != metric_filter.lower(): continue
            if noise_mode_filter and noise_mode != noise_mode_filter.lower(): continue
            if plevel_filter is not None and (not math.isfinite(plevel) or abs(plevel - plevel_filter) > 1e-12): continue
            rows.append(dict(metric=metric,k=k,L=L,n=n,noise=noise,plevel=plevel,noise_mode=noise_mode,
                             repeat=repeat,ARI=ARI,coherence_measured=coh))
    # filtre global : enlever K=1 (silhouette/indices pas comparables)
    rows = [r for r in rows if r["k"] != 1]
    return rows

def agg_mean(rows, keys, value_key):
    acc = defaultdict(list)
    for r in rows:
        v = r.get(value_key, float("nan"))
        if v == v:  # non-NaN
            acc[tuple(r[k] for k in keys)].append(v)
    return {k: float(np.mean(vs)) for k, vs in acc.items() if vs}

def sort_unique(vals): return sorted(set(vals))
def marker(i): return ["o","s","D","^","v",">","<","P","X","*"][i % 10]

def main():
    ap = argparse.ArgumentParser(description="Trace le dashboard ARI à partir de results_for_plot.csv")
    ap.add_argument("--results", default="results_for_plot.csv")
    ap.add_argument("--out", default="dashboard.png")
    ap.add_argument("--noise_ref", type=float, default=0.50)
    ap.add_argument("--noise_low", type=float, default=0.10)
    ap.add_argument("--noise_high", type=float, default=0.75)
    ap.add_argument("--metric_filter", default=None)         # ex: wmfd
    ap.add_argument("--noise_mode_filter", default=None)     # ex: drop / swap
    ap.add_argument("--plevel_filter", type=float, default=None) # ex: 0.30
    ap.add_argument("--use_measured_coherence", action="store_true")
    args = ap.parse_args()

    rows = load_rows(args.results,
                     metric_filter=args.metric_filter,
                     noise_mode_filter=args.noise_mode_filter,
                     plevel_filter=args.plevel_filter)

    if not rows:
        raise SystemExit("Aucune ligne exploitable. Vérifie --metric_filter/--noise_mode_filter/--plevel_filter et le CSV.")

    fig = plt.figure(figsize=(14,10), dpi=140)
    gs = fig.add_gridspec(2, 3, wspace=0.25, hspace=0.35)

    # (a) ARI vs noise (par K)
    ax = fig.add_subplot(gs[0,0]); ax.set_title("(a)", loc="left")
    ax.set_xlabel("noise (%)"); ax.set_ylabel("ARI")
    m = agg_mean(rows, keys=["k","noise"], value_key="ARI")
    Ks = sort_unique([k for (k,_n) in m.keys()])
    for i,k in enumerate(Ks):
        xs, ys = [], []
        for (kk, nv) in sorted(m.keys(), key=lambda t:(t[0], t[1])):
            if kk != k: continue
            xs.append(int(round(100*nv))); ys.append(m[(kk, nv)])
        if xs: ax.plot(xs, ys, marker=marker(i), label=f"K={k}")
    ax.set_ylim(0,1); ax.grid(True, alpha=.3)
    if Ks: ax.legend(frameon=False, fontsize=8)

    # Sous-ensemble à noise_ref
    rows_ref = [r for r in rows if math.isfinite(r["noise"]) and abs(r["noise"] - args.noise_ref) < 1e-12]

    # (b) ARI vs L @ noise_ref
    ax = fig.add_subplot(gs[0,1]); ax.set_title("(b)", loc="left")
    ax.set_xlabel("leaves (L)"); ax.set_ylabel("ARI")
    m = agg_mean(rows_ref, keys=["L"], value_key="ARI")
    xs = sorted(m.keys())
    if xs:
        ax.plot([x[0] for x in xs], [m[x] for x in xs], marker="o")
    ax.set_ylim(0,1); ax.grid(True, alpha=.3)

    # (c) ARI vs n @ noise_ref
    ax = fig.add_subplot(gs[0,2]); ax.set_title("(c)", loc="left")
    ax.set_xlabel("trees per group (n)"); ax.set_ylabel("ARI")
    m = agg_mean(rows_ref, keys=["n"], value_key="ARI")
    xs = sorted(m.keys())
    if xs:
        ax.plot([x[0] for x in xs], [m[x] for x in xs], marker="o")
    ax.set_ylim(0,1); ax.grid(True, alpha=.3)

    # (d) ARI vs plevel (ou cohérence mesurée)
    ax = fig.add_subplot(gs[1,0]); ax.set_title("(d)", loc="left")
    ax.set_ylabel("ARI"); ax.set_ylim(0,1); ax.grid(True, alpha=.3)
    if args.use_measured_coherence and any(math.isfinite(r["coherence_measured"]) for r in rows):
        ax.set_xlabel("coherence measured (%)")
        m = agg_mean(rows, keys=["coherence_measured"], value_key="ARI")
        xs = sorted([x for (x,) in m.keys() if x==x])
        if xs:
            ax.plot([int(round(100*x)) for x in xs], [m[(x,)] for x in xs], marker="o")
    else:
        ax.set_xlabel("plevel (%)")
        m = agg_mean(rows, keys=["plevel"], value_key="ARI")
        xs = sorted([x for (x,) in m.keys() if x==x])
        if xs:
            ax.plot([int(round(100*x)) for x in xs], [m[(x,)] for x in xs], marker="o")

    # (e) ARI vs K @ noise_ref
    ax = fig.add_subplot(gs[1,1]); ax.set_title("(e)", loc="left")
    ax.set_xlabel("K"); ax.set_ylabel("ARI")
    m = agg_mean(rows_ref, keys=["k"], value_key="ARI")
    xs = sorted(m.keys())
    if xs:
        ax.plot([x[0] for x in xs], [m[x] for x in xs], marker="o")
    ax.set_ylim(0,1); ax.grid(True, alpha=.3)

    # (f) ΔARI (low − high) par K
    ax = fig.add_subplot(gs[1,2]); ax.set_title("(f)", loc="left")
    ax.set_xlabel("K"); ax.set_ylabel("ΔARI")
    rows_low  = [r for r in rows if math.isfinite(r["noise"]) and abs(r["noise"] - args.noise_low)  < 1e-12]
    rows_high = [r for r in rows if math.isfinite(r["noise"]) and abs(r["noise"] - args.noise_high) < 1e-12]
    m_low  = agg_mean(rows_low,  keys=["k"], value_key="ARI")
    m_high = agg_mean(rows_high, keys=["k"], value_key="ARI")
    Ks_all = sorted(set([k for (k,) in m_low.keys()] + [k for (k,) in m_high.keys()]))
    if Ks_all:
        ys = []
        for k in Ks_all:
            a = m_low.get((k,), float("nan"))
            b = m_high.get((k,), float("nan"))
            ys.append((a - b) if (a==a and b==b) else float("nan"))
        ax.plot(Ks_all, ys, marker="o")
    ax.grid(True, alpha=.3)

    plt.savefig(args.out, bbox_inches="tight")
    print(f"[OK] Wrote {args.out}")

if __name__ == "__main__":
    main()