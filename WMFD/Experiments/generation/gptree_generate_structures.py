#!/usr/bin/env python3
# -- coding: utf-8 --
"""
gptree_generate_structures.py
-----------------------------
Génère en mémoire des "runs" d'arbres phylo-like (ETE3) selon des grilles de paramètres.
- Chaque run : K prototypes ; pour chaque prototype, n_per_group arbres perturbés.
- Perturbations = NNI + jitter longueurs ; Bruit = NNI additionnels + renommage de feuilles.
- Retourne des dataclasses (RunSpec) avec : trees_ete, trees_nested (optionnel), true_labels.

Dépendances : pip install ete3
"""

from __future__ import annotations
import random
from dataclasses import dataclass
from typing import List, Literal, Optional, Dict, Any, Tuple
from ete3 import Tree

# ======================= Structures de données =======================

@dataclass(frozen=True)
class Node:
    """Arbre immuable 100% Python (pas d'ETE3 requis côté consommateur)."""
    name: Optional[str]
    dist: float
    children: List["Node"]

@dataclass
class RunSpec:
    """Un 'run' avec méta + arbres (formats au choix)."""
    run_name: str
    K: int
    L: int
    n_per_group: int
    plevel: float
    noise_pct: float
    rep: int
    trees_ete: Optional[List[Tree]]           # présent si return_format inclut 'ete'
    trees_nested: Optional[List[Node]]        # présent si return_format inclut 'nested'
    true_labels: Optional[List[int]] = None   # étiquette vraie de cluster pour chaque arbre

ReturnFormat = Literal["ete", "nested", "both"]

# ======================= Génération des arbres =======================

def _make_base_tree(num_leaves: int, seed: Optional[int]=None) -> Tree:
    """
    Arbre binaire simple avec feuilles L1..Ln, longueurs > 0.
    """
    if seed is not None:
        random.seed(seed)

    t = Tree()
    r = t.add_child()
    a = r.add_child(name="L1"); a.dist = 0.1 + 0.2*random.random()
    b = r.add_child(name="L2"); b.dist = 0.1 + 0.2*random.random()

    for i in range(3, num_leaves + 1):
        internals = [n for n in t.traverse() if not n.is_leaf()]
        parent = random.choice(internals)
        new_internal = parent.add_child()
        if parent.children:
            ch = random.choice(parent.children)
            ch.detach()
            new_internal.add_child(ch)
        leaf = new_internal.add_child(name=f"L{i}")
        for node in (new_internal, leaf):
            node.dist = 0.1 + 0.2*random.random()
    return t

def _random_nni_moves(t: Tree, n_moves: int):
    """
    Petits NNI pour varier la topologie (sans supprimer de feuilles).
    """
    for _ in range(max(0, n_moves)):
        internals = [n for n in t.traverse()
                     if not n.is_leaf() and n.up is not None and len(n.children) >= 2]
        if not internals:
            return
        x = random.choice(internals)
        if x.up is None or len(x.up.children) < 2:
            continue
        sibs = [c for c in x.up.children if c is not x]
        if not sibs:
            continue
        a = random.choice(x.children)
        s = random.choice(sibs)
        b = random.choice(s.children) if s.children else s
        if a is None or b is None:
            continue
        pa, pb = a.up, b.up
        a.detach(); b.detach()
        pa.add_child(b); pb.add_child(a)

def _rescale_branch_lengths(t: Tree, factor: float = 1.0, jitter: float = 0.1):
    """
    Rescale + jitter contrôlé ; garantit dist >= 0.01.
    """
    for n in t.traverse():
        d = (n.dist if n.dist is not None else 0.1)
        d = d * factor * (1.0 + jitter*(random.random()-0.5))
        n.dist = max(0.01, d)

def _apply_noise_to_leafnames(t: Tree, noise_pct: float, seed: Optional[int]=None):
    """
    Bruit sur l'ensemble de taxons : remplace ~p% des feuilles par des noms uniques.
    Réduit l'intersection des noms entre arbres -> utile pour faire baisser ARI.
    """
    if noise_pct <= 0:
        return
    if seed is not None:
        random.seed(seed + 1_234_567)
    leaves = t.get_leaves()
    k = max(0, int(round(len(leaves) * (noise_pct/100.0))))
    idxs = random.sample(range(len(leaves)), k) if k>0 else []
    for i in idxs:
        leaves[i].name = f"NOISE_{i}_{random.randint(0, 10**9)}"

def _calc_internal_edges(t: Tree) -> int:
    return max(1, sum(1 for n in t.traverse() if not n.is_leaf() and n.up is not None))

def _make_tree_from_proto(proto: Tree, plevel: float, noise_pct: float=0.0, seed: Optional[int]=None) -> Tree:
    """
    Copie + perturbations contrôlées :
      - NNI (intensité ~ (1-plevel) + bruit)
      - Rescale longueurs (jitter ~ (1-plevel) + bruit)
      - Bruit sur noms de feuilles (proportion=p)
    Aucun changement du nombre de feuilles.
    """
    if seed is not None:
        random.seed(seed)
    t = proto.copy(method="deepcopy")

    internal_edges = _calc_internal_edges(t)
    base_moves  = int(round((1.0 - plevel) * internal_edges))
    extra_moves = int(round((noise_pct/100.0) * internal_edges))
    total_moves = base_moves + extra_moves
    _random_nni_moves(t, total_moves)

    jitter = 0.2*(1.0 - 0.5*plevel) + 0.3*(noise_pct/100.0)
    _rescale_branch_lengths(t, factor=1.0, jitter=jitter)

    _apply_noise_to_leafnames(t, noise_pct=noise_pct, seed=seed)
    return t

# ======================= Conversions =======================

def _ete_to_nested(n: Tree) -> Node:
    """
    Convertit un nœud ETE3 (TreeNode) en structure Node (récursif).
    """
    children = [ _ete_to_nested(c) for c in n.children ]
    name = n.name if (n.name not in (None, "", "NoName")) else None
    dist = float(n.dist) if n.dist is not None else 0.01
    return Node(name=name, dist=dist, children=children)

# ======================= Validations =======================

def validate_tree_ete(t: Tree) -> Dict[str, Any]:
    """
    Retourne un petit rapport de validation pour un arbre ETE3.
    """
    leaves = t.get_leaves()
    zero = sum(1 for n in t.traverse() if (n.dist is None or n.dist == 0))
    neg  = sum(1 for n in t.traverse() if (n.dist is not None and n.dist < 0))
    unlabeled = sum(1 for lf in leaves if not lf.name or lf.name.strip() == "")
    return {
        "num_leaves": len(leaves),
        "zero_lengths": zero,
        "negative_lengths": neg,
        "unlabeled_leaves": unlabeled,
        "ok": (zero == 0 and neg == 0 and unlabeled == 0 and len(leaves) >= 2),
    }

def validate_run(run: RunSpec) -> Dict[str, Any]:
    """
    Validation d’un run (agrégé).
    """
    reports = []
    if run.trees_ete:
        for t in run.trees_ete:
            reports.append(validate_tree_ete(t))
    ok = all(r["ok"] for r in reports) if reports else True
    return {"run_name": run.run_name, "trees": len(reports), "ok": ok, "details": reports[:5]}

# ======================= Générateur principal =======================

def generate_runs(
    Ks: List[int] = [1, 2, 3, 4],
    Ls: List[int] = [10, 20, 30, 40, 50, 60],
    ns: List[int] = [8, 16],
    plevels: List[float] = [0.30, 0.50, 0.70],
    noises: List[int] = [0, 25, 50, 75],
    reps: List[int] = [0, 1, 2],
    return_format: ReturnFormat = "both",
) -> List[RunSpec]:
    """
    Génère une liste de RunSpec. Aucun fichier n’est écrit ; tout est en mémoire.
    """
    all_runs: List[RunSpec] = []

    for K in Ks:
        for L in Ls:
            for n_per_group in ns:
                for plevel in plevels:
                    for noise in noises:
                        for rep in reps:
                            run_name = f"run_k{K}L{L}_n{n_per_group}_plev{int(plevel*100)}_p{int(noise)}rep{rep}"
                            trees_ete: Optional[List[Tree]] = [] if return_format in ("ete","both") else None
                            trees_nested: Optional[List[Node]] = [] if return_format in ("nested","both") else None
                            labels_true: List[int] = []

                            # K prototypes; chaque proto => n_per_group arbres perturbés
                            for g in range(K):
                                # NB: seed rend réplicable par (rep, g, L)
                                proto = _make_base_tree(L, seed=(rep*10000 + g*1000 + L))
                                for i in range(n_per_group):
                                    t = _make_tree_from_proto(
                                        proto,
                                        plevel=plevel,
                                        noise_pct=float(noise),
                                        seed=(rep*10000 + g*1000 + i)
                                    )
                                    if trees_ete is not None:
                                        trees_ete.append(t)
                                    if trees_nested is not None:
                                        root = t
                                        if len(root.children) == 1:
                                            root = root.children[0]
                                        trees_nested.append(_ete_to_nested(root))
                                    labels_true.append(f"{g}_{i}")   # g = cluster, i = index arbre dans ce cluster

                            all_runs.append(RunSpec(
                                run_name=run_name,
                                K=K, L=L, n_per_group=n_per_group,
                                plevel=plevel, noise_pct=float(noise),
                                rep=rep,
                                trees_ete=trees_ete,
                                trees_nested=trees_nested,
                                true_labels=labels_true,
                            ))
    return all_runs

# ======================= Démo / Test local =======================

if __name__ == "__main__":
    # Petit test autonome : génère quelques runs, imprime des aperçus (AUCUN fichier).
    runs = generate_runs(
        Ks=[1,2,3,4],                # 4 loops "K"
        Ls=[10, 20, 30, 40, 50, 70],             # 4 loops "L"
        ns=[8,16],                  # 4 loops "n"
        plevels=[0.3, 0.7],      # cohésion intra
        noises=[0, 25, 50, 75],          # 4 loops "p" (bruit)
        reps=[0,1],                # réplications
        return_format="both"     # "ete" | "nested" | "both"
    )
    print("Total runs:", len(runs), flush=True)

    # Validation rapide + aperçu texte de quelques arbres
    for r in runs:
        v = validate_run(r)
        print(f"- {r.run_name} | trees={len(r.trees_ete or [])} | valid={v['ok']}", flush=True)
        for t in (r.trees_ete or [])[:2]:
            newick = t.write(format=1)
            print("  newick:", (newick[:120] + "...") if len(newick) > 120 else newick, flush=True)

    # Exemple : accès aux labels vrais et à la structure imbriquée
    example = runs[0]
    print("\nExample run:", example.run_name, flush=True)
    print("true_labels (counts):", {lab: example.true_labels.count(lab) for lab in set(example.true_labels or [])}, flush=True)
    if example.trees_nested:
        print("nested sample:", example.trees_nested[0], flush=True)