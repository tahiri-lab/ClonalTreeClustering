# gptree_generate_structures.py
# Génère des runs "en mémoire" (AUCUN fichier).
# Retourne des dataclasses avec les arbres au format ETE3 et/ou structure imbriquée Python.
#
# Dépendances : pip install ete3

from __future__ import annotations
import random
from dataclasses import dataclass
from typing import List, Literal, Optional, Dict, Any
from ete3 import Tree

# ----------------------------- Structures de données -----------------------------

@dataclass
class Node:
    """Structure d'arbre immuable 100% Python (pas d'ETE3 requis côté consommateur)."""
    name: Optional[str]
    dist: float
    children: List["Node"]

@dataclass
class RunSpec:
    """Un 'run' avec méta + arbres (formats au choix, voir generate_runs)."""
    run_name: str
    K: int
    L: int
    n_per_group: int
    plevel: float
    noise_pct: float
    rep: int
    trees_ete: Optional[List[Tree]]           # présent si format demandé inclut 'ete'
    trees_nested: Optional[List[Node]]        # présent si format demandé inclut 'nested'
    # On peut ajouter d'autres vues (ex: newick str) si besoin sans écrire de fichiers.

# ----------------------------- Génération des arbres -----------------------------

def _make_base_tree(num_leaves: int, seed: Optional[int]=None) -> Tree:
    """Arbre binaire simple avec feuilles L1..Ln, longueurs > 0."""
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
    """Petits NNI pour varier la topologie (sans longueurs nulles)."""
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
    """Garantit dist >= 0.01 pour éviter les longueurs nulles."""
    for n in t.traverse():
        d = (n.dist if n.dist is not None else 0.1)
        d = d * factor * (1.0 + jitter*(random.random()-0.5))
        n.dist = max(0.01, d)

def _make_tree_from_proto(proto: Tree, plevel: float, seed: Optional[int]=None) -> Tree:
    """Copie + perturbations contrôlées (NNI + rescale). Aucune suppression de feuilles."""
    if seed is not None:
        random.seed(seed)
    t = proto.copy(method="deepcopy")
    internal_edges = max(1, sum(1 for n in t.traverse() if not n.is_leaf() and n.up is not None))
    n_moves = int(round((1.0 - plevel) * internal_edges))
    _random_nni_moves(t, n_moves)
    _rescale_branch_lengths(t, factor=1.0, jitter=0.2*(1.0 - 0.5*plevel))
    return t

# ----------------------------- Conversions ---------------------------------------

def _ete_to_nested(n: Tree) -> Node:
    """Convertit un nœud ETE3 (TreeNode) en structure Node (récursif)."""
    children = [ _ete_to_nested(c) for c in n.children ]
    name = n.name if (n.name not in (None, "", "NoName")) else None
    dist = float(n.dist) if n.dist is not None else 0.01
    return Node(name=name, dist=dist, children=children)

# ----------------------------- Validations (in-memory) ---------------------------

def validate_tree_ete(t: Tree) -> Dict[str, Any]:
    """Retourne un petit rapport de validation pour un arbre ETE3."""
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
    """Validation d’un run (agrégé)."""
    reports = []
    if run.trees_ete:
        for t in run.trees_ete:
            reports.append(validate_tree_ete(t))
    ok = all(r["ok"] for r in reports) if reports else True
    return {"run_name": run.run_name, "trees": len(reports), "ok": ok, "details": reports[:5]}  # aperçus

# ----------------------------- Générateur principal ------------------------------

ReturnFormat = Literal["ete", "nested", "both"]

def generate_runs(
    Ks: List[int] = [1, 2, 3, 4],
    Ls: List[int] = [10, 20, 30, 40, 50, 60],
    ns: List[int] = [8, 16],
    plevels: List[float] = [0.30, 0.50, 0.70],
    noises: List[int] = [0, 25, 50, 75],        # gardé pour compat signature; pas utilisé pour "drop"
    reps: List[int] = [0, 1, 2],
    return_format: ReturnFormat = "both",
) -> List[RunSpec]:
    """
    Génère tous les runs en mémoire et renvoie une liste de RunSpec.
    Aucun fichier n'est écrit.
    - return_format: 'ete' (objets ETE3), 'nested' (Node), ou 'both'.
    """
    all_runs: List[RunSpec] = []

    for K in Ks:
        for L in Ls:
            for n_per_group in ns:
                for plevel in plevels:
                    for noise in noises:
                        for rep in reps:
                            run_name = f"run_k{K}L{L}_n{n_per_group}_plev{int(plevel*100)}_p{int(noise)}rep{rep}"
                            trees_ete: Optional[List[Tree]] = [] if return_format in ("ete", "both") else None
                            trees_nested: Optional[List[Node]] = [] if return_format in ("nested", "both") else None

                            # K prototypes; chaque proto => n_per_group arbres perturbés
                            for g in range(K):
                                proto = _make_base_tree(L, seed=(rep*10000 + g*1000 + L))
                                for i in range(n_per_group):
                                    t = _make_tree_from_proto(proto, plevel=plevel, seed=(rep*10000 + g*1000 + i))
                                    if trees_ete is not None:
                                        trees_ete.append(t)
                                    if trees_nested is not None:
                                        # convertir à partir de la racine implicite de ETE3
                                        # Attention: ETE3 a une racine implicite; on descend si un seul enfant non-feuille
                                        root = t
                                        if len(root.children) == 1:
                                            root = root.children[0]
                                        trees_nested.append(_ete_to_nested(root))

                            all_runs.append(RunSpec(
                                run_name=run_name,
                                K=K, L=L, n_per_group=n_per_group,
                                plevel=plevel, noise_pct=float(noise),
                                rep=rep,
                                trees_ete=trees_ete,
                                trees_nested=trees_nested
                            ))
    return all_runs

# ----------------------------- Exemple d’utilisation -----------------------------

# gptree_generate_structures.py

# ... tout ton code des classes et fonctions ...

if __name__ == "__main__":
    
    runs = generate_runs(
        Ks=[1,2,3,4],
        Ls=[10, 20, 30, 40, 50, 60],
        ns=[8,16],
        plevels=[0.3,0.5,0.7],
        noises=[0,25,50,75],
        reps=[0,1,2],
        return_format="both"
    )
    print("Total runs:", len(runs))

    
    for r in runs:
        print(r.run_name, ":", len(r.trees_ete), "arbres")
        for t in r.trees_ete[:4]:
            s = t.write(format=1)
            print("  -", (s[:100]+"...") if len(s)>100 else s)