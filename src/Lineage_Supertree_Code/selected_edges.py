from ete3 import Tree
from collections import defaultdict, Counter
from math import ceil
from cluster import cluster_finder
from average_Jaccard_distance import compute_average_jaccard_distance
from frequent_arc import extract_frequent_arcs
from Reference_tab import build_reference_matrix_from_nonfrequent_arcs, get_direct_arcs
from group_edges import group_arcs_by_child

def build_global_arc_map(trees):
    arc_map = defaultdict(list)
    for i, newick in enumerate(trees):
        tree = Tree(newick, format=8)
        arcs = get_direct_arcs(tree)
        for arc in arcs:
            arc_map[arc].append(i)
    return arc_map

def select_arcs_by_criteria(grouped_arcs, cluster_indices, global_arc_map):
    final_selected = set()
    remaining_unresolved = set()
    threshold = ceil(len(cluster_indices) / 2)

    for child, arc_counts in grouped_arcs.items():
        # 1. Garder les arcs qui passent le seuil de fréquence
        valid_arcs = [(arc, count) for arc, count in arc_counts.items() if count >= threshold-1]
        if not valid_arcs:
            continue

        # 2. Sélectionner ceux avec le support max dans le cluster
        max_count = max(count for arc, count in valid_arcs)
        candidates = [arc for arc, count in valid_arcs if count == max_count]

        if len(candidates) == 1:
            final_selected.add(candidates[0])
        else:
            # 3. Appliquer critère 2 : le support externe minimal
            external_supports = {
                arc: len(set(global_arc_map[arc]) - set(cluster_indices)) for arc in candidates
            }
            min_support = min(external_supports.values())
            best_arcs = [arc for arc in candidates if external_supports[arc] == min_support]

            if len(best_arcs) == 1:
                final_selected.add(best_arcs[0])
            else:
                # Cas d’égalité non résolu
                remaining_unresolved.update(best_arcs)

    return final_selected, remaining_unresolved

def build_supertrees_by_cluster(trees, clusters, frequent_arcs, global_arc_map):
    supertrees = []
    unselected_arcs = []

    for cluster in clusters:
        cluster_trees = [trees[i] for i in cluster]
        grouped_arcs = group_arcs_by_child(cluster_trees, frequent_arcs)
        selected, not_selected = select_arcs_by_criteria(grouped_arcs, cluster, global_arc_map)

        full_tree_arcs = frequent_arcs.union(selected)
        supertrees.append(full_tree_arcs)
        unselected_arcs.append(not_selected)

    return supertrees, unselected_arcs

# MAIN
if __name__ == "__main__":
    trees = [
        "(((12)2)9,((6,1)7,(3,5,13)10)4,14,(11)15,(16)8)N;",
        "(((2)9)12,((1,4)7,(10,5,13)6)3,14,(15)11,(16)8)N;",
        "(((12)9)2,((13,5)1,(3,4,6)7)10,14,(11)15,(16)8)N;",
        "(((12)2)9,((3,1)7,(13,6,5)4)10,14,(11)15,(8)16)N;",
        "(((9)2)12,((3,1)7,(13,4,6)5)10,14,(11)15,(8)16)N;"
    ]

    threshold = compute_average_jaccard_distance(trees, verbose=False)
    frequent_arcs = extract_frequent_arcs(trees, threshold, verbose=False)
    global_arc_map = build_global_arc_map(trees)
    ref_matrix = build_reference_matrix_from_nonfrequent_arcs(trees, frequent_arcs, verbose=False)
    clusters = cluster_finder(ref_matrix)

    supertrees, unselected = build_supertrees_by_cluster(trees, clusters, frequent_arcs, global_arc_map)

    for i, arcs in enumerate(supertrees):
        print(f"\n✅ Cluster {i+1} – {len(arcs)} selected arcs:")
        for arc in sorted(arcs):
            print(f"  {arc}")

    for i, arcs in enumerate(unselected):
        print(f"\n❌ Cluster {i+1} – {len(arcs)} unresolved arcs:")
        for arc in sorted(arcs):
            print(f"  {arc}")