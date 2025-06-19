import re

# Expression régulière pour capturer name@weight:branch dans les lignes Newick
pattern = re.compile(r'([A-Za-z0-9_]+)@([0-9]+):([0-9]+)')

# Ouverture du fichier et extraction des 5 premiers arbres
with open('weighted_newicks_60.txt', 'r') as f:
    lines = [l.strip() for l in f if l.strip()]
first5 = lines[:5]

# Structures pour stocker les informations par arbre
tree_ids = []
tree_adjs = {}      # dictionnaire {arbre: liste d'adjacence}
tree_weights = {}   # dictionnaire {arbre: {noeud: poids}}
tree_branch = {}    # dictionnaire {arbre: {noeud: longueur_de_branche}}
tree_outdeg = {}    # dictionnaire {arbre: {noeud: degré de sortie}}
tree_leaves = {}    # dictionnaire {arbre: ensemble des feuilles}
tree_heights = {}   # dictionnaire {arbre: {noeud: hauteur}}

for line in first5:
    # Identifiant de l'arbre (ex: "60_1", "60_2", ...)
    tree_id = line.split(':', 1)[0].strip()
    tree_ids.append(tree_id)
    
    # Extraction des paires parent->enfants depuis le JSON pré-traité (adjacency_lists.json)
    # (Chaque arbre possède "naive" comme racine)
    # Ici, on suppose avoir un dictionnaire 'processed' pré-chargé contenant les listes d'adjacence
    adj = processed[tree_id]["adjacency_list"].copy()
    if "ROOT" in adj:
        adj.pop("ROOT")   # Retirer le noeud racine artificiel "ROOT" si présent
    tree_adjs[tree_id] = adj
    
    # Extraction des informations de poids et longueur de branche via regex
    matches = pattern.findall(line)
    weight_map = {name: int(wt) for name, wt, br in matches}
    branch_map = {name: int(br) for name, wt, br in matches if name != "naive"}
    tree_weights[tree_id] = weight_map
    tree_branch[tree_id] = branch_map
    
    # Calcul des degrés de sortie (nombre d'enfants) pour chaque nœud à partir de la liste d'adjacence
    outdeg_map = {}
    for parent, children in adj.items():
        outdeg_map[parent] = len(children)
        for child in children:
            outdeg_map.setdefault(child, 0)
    tree_outdeg[tree_id] = outdeg_map
    
    # Identification des feuilles (nœuds sans enfants)
    leaves = {node for node, deg in outdeg_map.items() if deg == 0}
    if "naive" in leaves:
        leaves.remove("naive")  # Exclure la racine "naive" des feuilles
    tree_leaves[tree_id] = leaves
    
    # Calcul de la hauteur de chaque nœud (distance en arêtes jusqu'à la racine "naive")
    height_map = {}
    from collections import deque
    root = "naive"
    queue = deque([(root, 1)])  # hauteur de la racine = 1 (par définition de l'indice de hauteur)
    visited = {root}
    while queue:
        node, h = queue.popleft()
        height_map[node] = h
        for child in adj.get(node, []):
            if child not in visited:
                visited.add(child)
                queue.append((child, h + 1))
    tree_heights[tree_id] = height_map

# Vérification du nombre de nœuds et de feuilles extraits pour chaque arbre
for tid in tree_ids:
    print(f"{tid}: {len(tree_weights[tid])} nœuds au total, dont {len(tree_leaves[tid])} feuilles.")
