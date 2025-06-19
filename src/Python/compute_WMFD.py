# Définition des poids pour les 5 features (branch length, height, weight, out-degree, topology)
lambda1, lambda2, lambda3, lambda4, lambda5 = 0.30, 0.20, 0.25, 0.15, 0.10

# Initialisation d'une matrice de distances 5x5 (remplie symétriquement)
n = len(tree_ids)
dist_matrix = [[0.0] * n for _ in range(n)]

# Pré-calcul des ensembles d'arêtes pour chaque arbre (pour accélérer le calcul de HD)
edges_per_tree = {}
for idx, tid in enumerate(tree_ids):
    edge_set = set()
    for parent, children in tree_adjs[tid].items():
        for child in children:
            edge_set.add((parent, child))
    edges_per_tree[tid] = edge_set

# Calcul de la distance pour chaque paire (i, j)
for i in range(n):
    T1 = tree_ids[i]
    for j in range(i, n):
        T2 = tree_ids[j]
        if i == j:
            dist_matrix[i][j] = 0.0
        else:
            # 1. Nœuds communs et distincts
            common_leaves = tree_leaves[T1] & tree_leaves[T2]
            union_leaves = tree_leaves[T1] | tree_leaves[T2]
            CN = len(common_leaves)  # nombre de feuilles communes
            TN = len(union_leaves)    # nombre total de feuilles uniques
            P = 1 - (CN / TN) if TN > 0 else 0  # indice de pénalité
            
            # Sommes des différences pour features
            common_sum_BL = common_sum_H = common_sum_W = common_sum_D = 0.0
            uncommon_sum_BL = uncommon_sum_H = uncommon_sum_W = uncommon_sum_D = 0.0
            
            # 2 & 3. Calcul des différences pour chaque nœud de l'union des feuilles
            for node in union_leaves:
                # Valeurs des features dans chaque arbre (0 si le nœud est absent)
                bl1 = tree_branch[T1].get(node, 0)
                bl2 = tree_branch[T2].get(node, 0)
                h1  = tree_heights[T1].get(node, 0)
                h2  = tree_heights[T2].get(node, 0)
                w1  = tree_weights[T1].get(node, 0)
                w2  = tree_weights[T2].get(node, 0)
                d1  = tree_outdeg[T1].get(node, 0)
                d2  = tree_outdeg[T2].get(node, 0)
                # Différences absolues pour chaque feature
                diff_BL = abs(bl1 - bl2)
                diff_H  = abs(h1 - h2)
                diff_W  = abs(w1 - w2)
                diff_D  = abs(d1 - d2)
                if node in common_leaves:
                    # Nœud commun aux deux arbres
                    common_sum_BL += diff_BL
                    common_sum_H  += diff_H
                    common_sum_W  += diff_W
                    common_sum_D  += diff_D
                else:
                    # Nœud distinct (présent dans un arbre seulement)
                    uncommon_sum_BL += diff_BL
                    uncommon_sum_H  += diff_H
                    uncommon_sum_W  += diff_W
                    uncommon_sum_D  += diff_D
            
            # Moyennes des différences (par nœud) pour chaque feature
            BL_common = common_sum_BL / CN if CN > 0 else 0.0
            H_common  = common_sum_H  / CN if CN > 0 else 0.0
            W_common  = common_sum_W  / CN if CN > 0 else 0.0
            D_common  = common_sum_D  / CN if CN > 0 else 0.0
            nb_uncommon = TN - CN  # nombre de nœuds distincts au total
            BL_uncommon = uncommon_sum_BL / nb_uncommon if nb_uncommon > 0 else 0.0
            H_uncommon  = uncommon_sum_H  / nb_uncommon if nb_uncommon > 0 else 0.0
            W_uncommon  = uncommon_sum_W  / nb_uncommon if nb_uncommon > 0 else 0.0
            D_uncommon  = uncommon_sum_D  / nb_uncommon if nb_uncommon > 0 else 0.0
            
            # 4. Combinaison pondérée des features (Weighted Node Differences)
            WND_common   = lambda1 * BL_common   + lambda2 * H_common   + lambda3 * W_common   + lambda4 * D_common
            WND_uncommon = lambda1 * BL_uncommon + lambda2 * H_uncommon + lambda3 * W_uncommon + lambda4 * D_uncommon
            
            # 5. Distance de Hamming pour la topologie (différences d'arêtes)
            edges_T1 = edges_per_tree[T1]
            edges_T2 = edges_per_tree[T2]
            union_edges = edges_T1 | edges_T2
            common_edges = edges_T1 & edges_T2
            HD = 1 - (len(common_edges) / len(union_edges)) if len(union_edges) > 0 else 0.0
            
            # 6. Calcul final de WMFD selon l'équation (14)
            WMFD_T1_T2 = P * WND_uncommon + WND_common + lambda5 * HD
            dist_matrix[i][j] = WMFD_T1_T2
            dist_matrix[j][i] = WMFD_T1_T2

# Affichage de la matrice des distances calculées (pour vérification)
for row in dist_matrix:
    print(row)
  import csv

with open('distances_WMFD.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["Tree"] + tree_ids)
    for i, tid in enumerate(tree_ids):
        row = [tid] + [f"{dist_matrix[i][j]:.6f}" for j in range(n)]
        writer.writerow(row)
print("✅ Fichier distances_WMFD.csv généré avec succès.")
