# -- coding: utf-8 --
"""
grid_experiments.py

Orchestrate the steps 1a 1c of the pipeline on a grid (n,k,L,p,noise,repeat),
then aggregates the results (ARI, objective, consistency) in a CSV.

✔ Loops in order:  n k L p  (and noise, repeat)
✔ --resume: skips runs already done (present in --out or metrics.csv)
✔ Generator: --gen_timeout_s (0 = no timeout), --gen_max_tries
✔ WMFD: configurable weights --wmfd_weights 0.30,0.20,0.25,0.15,0.10
✔ Intra-cluster consistency Jaccard: before and after noise
✔ Timings by step: sec_gen / sec_metric / sec_cluster

CSV output (minimal columns expected per plot):
timestamp,metric,k,L,n_per_group,noise,noise_mode,plevel,repeat,ARI,objective,coherence_measured,seed
+ bonus columns: run_dir, coherence_pre, sec_gen, sec_metric, sec_cluster
"""

import argparse, os, sys, subprocess, random, csv, pathlib, time
from datetime import datetime
from itertools import combinations
from typing import List

from ete3 import Tree

THIS_PY = sys.executable  # current interpreter (venv)


# ---------- FS / subprocess utils ----------

def ensure_dir(path: str):
    os.makedirs(path, exist_ok=True)

def run(cmd: List[str]):
    subprocess.run(cmd, check=True)

def parse_list(s, cast=float):
    return [cast(x.strip()) for x in s.split(",") if x.strip()]

def _timed(fn, *a):
    t0 = time.perf_counter()
    fn(*a)
    return time.perf_counter() - t0


# --------- noise & consistency ---------

def apply_label_noise_to_csv(in_csv, out_csv, noise_frac, seed):
    """Bruit 'swap' : échange ~noise_frac des feuilles (par paires) dans chaque arbre."""
    if noise_frac <= 0:
        with open(in_csv, "r", encoding="utf-8") as fsrc, open(out_csv, "w", encoding="utf-8", newline="") as fdst:
            fdst.write(fsrc.read())
        return
    rng = random.Random(seed)
    with open(in_csv, newline="", encoding="utf-8") as f_in, \
         open(out_csv, "w", newline="", encoding="utf-8") as f_out:
        r = csv.DictReader(f_in)
        w = csv.DictWriter(f_out, fieldnames=r.fieldnames)
        w.writeheader()
        for row in r:
            t = Tree(row["newick"].strip(), format=1)
            leaves = [lf for lf in t.iter_leaves()]
            L = len(leaves)
            swaps = max(0, int((noise_frac * L) // 2))
            if swaps > 0:
                idx = rng.sample(range(L), 2 * swaps)
                for a, b in zip(idx[::2], idx[1::2]):
                    leaves[a].name, leaves[b].name = leaves[b].name, leaves[a].name
            row["newick"] = t.write(format=1).strip()
            w.writerow(row)

def apply_leaf_drop_noise_to_csv(in_csv, out_csv, drop_frac, seed):
    """Bruit 'drop' : retire aléatoirement ~drop_frac des feuilles dans chaque arbre (≥2 feuilles conservées)."""
    if drop_frac <= 0:
        with open(in_csv, "r", encoding="utf-8") as fsrc, open(out_csv, "w", encoding="utf-8", newline="") as fdst:
            fdst.write(fsrc.read())
        return
    rng = random.Random(seed)
    with open(in_csv, newline="", encoding="utf-8") as f_in, \
         open(out_csv, "w", newline="", encoding="utf-8") as f_out:
        r = csv.DictReader(f_in)
        w = csv.DictWriter(f_out, fieldnames=r.fieldnames)
        w.writeheader()
        for row in r:
            t = Tree(row["newick"].strip(), format=1)
            leaf_names = [lf.name for lf in t.iter_leaves()]
            L = len(leaf_names)
            m = int(round(drop_frac * L))
            m = max(0, min(L - 2, m))  # keep at least 2 sheets
            if m > 0:
                to_drop = set(rng.sample(leaf_names, m))
                keep = [x for x in leaf_names if x not in to_drop]
                t.prune(keep, preserve_branch_length=True)
            row["newick"] = t.write(format=1).strip()
            w.writerow(row)

def measure_intra_cluster_coherence(csv_path):
    """Cohérence intra-cluster = Jaccard moyen des ensembles de feuilles sur toutes les paires d'un cluster."""
    by_cluster = {}
    with open(csv_path, newline="", encoding="utf-8") as f:
        r = csv.DictReader(f)
        for row in r:
            c = int(str(row["cluster"]).strip())
            leaves = set(Tree(row["newick"], format=1).get_leaf_names())
            by_cluster.setdefault(c, []).append(leaves)
    vals = []
    for sets in by_cluster.values():
        if len(sets) < 2:
            continue
        for a, b in combinations(sets, 2):
            u = len(a | b)
            vals.append((len(a & b) / u) if u else 0.0)
    return (sum(vals) / len(vals)) if vals else 0.0


# ---------- main ----------

def main():
    ap = argparse.ArgumentParser(description="Run (n,k,L,p,noise) grid and collect ARI for dashboard.")
    ap.add_argument("--workdir", default="C:\\ARI_Dashboard\\Essaie",
                    help="Dossier de travail pour créer les runs (run_*) et écrire les sorties.")
    ap.add_argument("--scripts_dir", default=None,
                    help="Dossier où se trouvent les scripts 1a/1b/1c. Par défaut = --workdir.")

    ap.add_argument("--k_list",   default="3",           help="ex: 2,3,4,5")
    ap.add_argument("--L_list",   default="20,40,60",    help="ex: 20,40,60")
    ap.add_argument("--n_list",   default="8,16,32,64",  help="nb arbres/cluster")
    ap.add_argument("--plevel_list", default="0.3,0.5,0.6", help="cohérence intra-cluster cible du générateur")
    ap.add_argument("--noise_list",  default="0.10,0.25,0.50,0.75", help="fraction de bruit (swap/drop)")
    ap.add_argument("--noise_mode", choices=["swap", "drop"], default="swap",
                    help="Type de bruit appliqué après génération.")

    ap.add_argument("--metric", choices=["wmfd", "rf"], default="wmfd")
    ap.add_argument("--wmfd_weights", default="0.30,0.20,0.25,0.15,0.10",
                    help="Poids L1..L5 pour WMFD, ex: 0.30,0.20,0.25,0.15,0.10")

    ap.add_argument("--repeats", type=int, default=1)
    ap.add_argument("--seed",    type=int, default=0)
    ap.add_argument("--resume", action="store_true",
                    help="Saute les runs déjà présents dans --out ou déjà complétés (metrics.csv).")

    ap.add_argument("--gen_timeout_s",   type=int, default=0,    help="Timeout par cluster pour générateur (0 = off).")
    ap.add_argument("--gen_max_tries",   type=int, default=2000, help="Max tentatives/ajout d'arbre (générateur).")

    ap.add_argument("--out", default="results_dashboard.csv", help="CSV d’agrégation des résultats.")
    args = ap.parse_args()

    workdir = os.path.abspath(args.workdir)
    scripts_dir = os.path.abspath(args.scripts_dir) if args.scripts_dir else workdir
    ensure_dir(workdir)

    Ks  = parse_list(args.k_list, int)
    Ls  = parse_list(args.L_list, int)
    Ns  = parse_list(args.n_list, int)
    PLs = parse_list(args.plevel_list, float)
    Ps  = parse_list(args.noise_list, float)
    L1, L2, L3, L4, L5 = [float(x) for x in args.wmfd_weights.split(",")]

    # scripts 
    path_gen   = os.path.join(scripts_dir, "gptree_cluster_refined.py")
    path_m_wmf = os.path.join(scripts_dir, "step1b_metric_wmfd.py")
    path_m_rf  = os.path.join(scripts_dir, "step1b_metric_rf.py")
    path_km    = os.path.join(scripts_dir, "step1c_clustering_kmedoids.py")

    # quick checks
    if not os.path.isfile(path_gen):
        print(f"[FATAL] Introuvable: {path_gen}"); sys.exit(2)
    if args.metric == "wmfd":
        if not os.path.isfile(path_m_wmf):
            print(f"[FATAL] Introuvable: {path_m_wmf}"); sys.exit(2)
    else:
        if not os.path.isfile(path_m_rf):
            print(f"[FATAL] Introuvable: {path_m_rf}"); sys.exit(2)
    if not os.path.isfile(path_km):
        print(f"[FATAL] Introuvable: {path_km}"); sys.exit(2)

    # CSV output
    out_csv = os.path.abspath(args.out)
    ensure_dir(str(pathlib.Path(out_csv).parent))
    first = not os.path.exists(out_csv)

    # recovery: tags already present
    done_tags = set()
    if args.resume and os.path.exists(out_csv):
        with open(out_csv, newline="", encoding="utf-8") as fprev:
            for r in csv.DictReader(fprev):
                try:
                    k0   = int(r["k"]); L0 = int(r["L"]); n0 = int(r["n_per_group"])
                    p0   = float(r["plevel"]); ptag = int(round(100 * p0))
                    nse0 = float(r["noise"]);  ntag = int(round(100 * nse0))
                    tag0 = f"k{k0}L{L0}_n{n0}_plev{ptag}_p{ntag}{r['noise_mode']}_rep{r['repeat']}"
                    done_tags.add(tag0)
                except Exception:
                    pass
        print(f"[RESUME] {len(done_tags)} runs déjà présents, ils seront ignorés.")

    with open(out_csv, "a", newline="", encoding="utf-8") as f_out:
        w = csv.writer(f_out)
        if first:
            w.writerow([
                "timestamp","metric","k","L","n_per_group","noise","noise_mode","plevel","repeat",
                "ARI","objective","coherence_measured","seed",
                "run_dir","coherence_pre","sec_gen","sec_metric","sec_cluster"
            ])

        base_seed = args.seed

        # ----- 4 LOOPS: (n, k, L, p) -----
        for n_per in Ns:                 # 1) n (trees/cluster)
            for k in Ks:                 # 2) k
                for L in Ls:             # 3) L (leaves)
                    for plevel in PLs:   # 4) p (overlap target)
                        for noise in Ps:
                            for rep in range(args.repeats):
                                tag = f"k{k}L{L}_n{n_per}_plev{int(round(plevel*100))}_p{int(round(noise*100))}{args.noise_mode}_rep{rep}"
                                run_dir = os.path.join(workdir, f"run_{tag}")
                                ensure_dir(run_dir)

                                # skip if (CSV) or if metrics.csv present
                                if args.resume and (tag in done_tags or os.path.exists(os.path.join(run_dir, "metrics.csv"))):
                                    print(f"[SKIP resume] {tag}")
                                    continue

                                seed = base_seed + (hash(("v2", n_per, k, L, plevel, noise, rep, args.noise_mode, args.metric)) & 0xffff)

                                trees_csv   = os.path.join(run_dir, "trees.csv")
                                noisy_csv   = os.path.join(run_dir, "trees_noisy.csv")
                                matrix_csv  = os.path.join(run_dir, "matrix.csv")
                                labels_csv  = os.path.join(run_dir, "labels.csv")
                                pairs_csv   = os.path.join(run_dir, "pairs.csv")
                                metrics_csv = os.path.join(run_dir, "metrics.csv")

                                print(f"\n=== (n={n_per}, k={k}, L={L}, p={plevel:.2f}) noise={noise:.2f} rep={rep} ===")

                                # 1a) generator
                                cmd_gen = [
                                    THIS_PY, path_gen,
                                    "--k", str(k), "--L", str(L), "--Ngen", str(n_per),
                                    "--plevel", str(plevel), "--seed", str(seed),
                                    "--timeout_s", str(args.gen_timeout_s),
                                    "--max_tries_per_tree", str(args.gen_max_tries),
                                    "--out", trees_csv
                                ]
                                print(f"[GEN] {tag}")
                                t_gen = _timed(run, cmd_gen)

                                # noise
                                if args.noise_mode == "swap":
                                    apply_label_noise_to_csv(trees_csv, noisy_csv, noise, seed)
                                else:
                                    apply_leaf_drop_noise_to_csv(trees_csv, noisy_csv, noise, seed)

                                # coherences
                                coh_pre  = measure_intra_cluster_coherence(trees_csv)
                                coh_post = measure_intra_cluster_coherence(noisy_csv)
                                print(f"[COH] target(p)={plevel:.2f}  pre={coh_pre:.3f}  post={coh_post:.3f}")

                                # 1b) metric
                                if args.metric == "wmfd":
                                    cmd_m = [
                                        THIS_PY, path_m_wmf,
                                        "--in_csv", noisy_csv,
                                        "--out_matrix", matrix_csv,
                                        "--out_pairs", pairs_csv,
                                        "--out_labels", labels_csv,
                                        "--L1", str(L1), "--L2", str(L2),
                                        "--L3", str(L3), "--L4", str(L4), "--L5", str(L5),
                                    ]
                                else:
                                    cmd_m = [
                                        THIS_PY, path_m_wmf,
                                        "--in_csv", noisy_csv,
                                        "--out_matrix", matrix_csv,
                                        "--out_pairs", pairs_csv,
                                        "--out_labels", labels_csv,
                                    ]
                                print(f"[METRIC-{args.metric.upper()}] {tag}")
                                t_metric = _timed(run, cmd_m)

                                # 1c) clustering
                                cmd_c = [
                                    THIS_PY, path_km,
                                    "--matrix", matrix_csv, "--k", str(k),
                                    "--labels", labels_csv, "--seed", str(seed),
                                    "--n_init", "50", "--max_iter", "200",
                                    "--out_pred", os.path.join(run_dir, "clusters_pred.csv"),
                                    "--out_medoids", os.path.join(run_dir, "medoids.csv"),
                                    "--out_kxn", os.path.join(run_dir, "dist_to_medoids.csv"),
                                    "--out_metrics", metrics_csv
                                ]
                                print(f"[KMEDOIDS] {tag}")
                                t_cluster = _timed(run, cmd_c)

                                # Read metrics.csv
                                with open(metrics_csv, newline="", encoding="utf-8") as fm:
                                    r = next(csv.DictReader(fm))
                                # ARI
                                ari = r.get("ARI")
                                if ari is None:
                                    ari = r.get("ari", "nan")
                                ari = float(ari)
                                # objective
                                obj = r.get("objective", "")
                                obj = float(obj) if obj not in ("", None, "nan") else None

                                # write row
                                w.writerow([
                                    datetime.now().isoformat(timespec="seconds"),
                                    args.metric, k, L, n_per, noise, args.noise_mode,
                                    plevel, rep, ari, obj, coh_post, seed,
                                    run_dir, f"{coh_pre:.6f}", f"{t_gen:.3f}", f"{t_metric:.3f}", f"{t_cluster:.3f}"
                                ])
                                f_out.flush()
                                print(f"[DONE] {tag}  ARI={ari:.6f}  obj={obj}  "
                                      f"coh_pre={coh_pre:.3f} coh_post={coh_post:.3f}  "
                                      f"t(gen/metric/cluster)={t_gen:.1f}/{t_metric:.1f}/{t_cluster:.1f}s")

if __name__ == "__main__":
    main()