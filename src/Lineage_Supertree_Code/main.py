from Newick_reader import *
from average_Jaccard_distance import *
from frequent_arc import *
from jaccard_distance import *
from Reference_tab import *
from Ancestral_distance import *
from cluster import *
from degree_of_membership import *
from doubted_edges import *
from group_edges import *
from probability_tab import *
from selected_edges import *

# --------------------------------------------------------
# Step 1: Load and validate input trees in Newick format
# --------------------------------------------------------

# Define the set of trees to process (can also be loaded from files)
trees = [
    "(((12)2)9,((1,4)7,(5,6,13)3)10,14,(11)15,(16)8)N;",
    "(((2)9)12,((3,1)7,(13,5,6)4)10,14,(15)11,(8)16)N;",
    "(((12)2)9,((13,4)7,(5,6,1)3)10,14,(15)11,(16)8)N;",
    "(((9)2)12,((1,3)7,(4,5,13)6)10,14,(11)15,(8)16)N;",
    "(((12)9)2,((6,1)7,(13,4,5)3)10,14,(11)15,(16)8)N;"
]

# Example of loading and validating a tree from a file
# newick_tree = read_newick(file_path)
# trees.append(newick_tree)

# ----------------------------------------------------------------
# Step 2: Compute the average Jaccard distance between all trees
#         This threshold is used to determine globally frequent arcs
# ----------------------------------------------------------------
threshold = compute_average_jaccard_distance(trees)

# --------------------------------------------------------
# Step 3: Extract arcs that are globally frequent
#         These arcs appear consistently across the tree set
# --------------------------------------------------------
frequents_arcs = extract_frequent_arcs(trees, threshold)

# ------------------------------------------------------------------
# Step 4: Build a reference matrix using only non-frequent arcs
#         This matrix captures how many arcs each pair of trees shares
# ------------------------------------------------------------------
reference_matrix = build_reference_matrix_from_nonfrequent_arcs(trees, frequents_arcs)

# --------------------------------------------------------------
# Step 5: Detect clusters of trees based on shared arc patterns
# --------------------------------------------------------------
clusters = cluster_finder(reference_matrix)
print("\n📊 Identified clusters of similar trees:", clusters)

# -----------------------------------------------------------------------
# Step 6: Build supertrees from each cluster using frequent & selected arcs
#         Also identify arcs that could not be reliably assigned
# -----------------------------------------------------------------------
global_arc_map = build_global_arc_map(trees)

# Returns: supertrees = list of arc sets (one per cluster),
#          unselected = list of unresolved arcs (one per cluster)
supertrees, unselected = build_supertrees_by_cluster(
    trees, clusters, frequents_arcs, global_arc_map
)

# -------------------------------------------------------
# Step 7: Print results – selected and unresolved arcs
# -------------------------------------------------------
for i, arcs in enumerate(supertrees):
    print(f"\n✅ Cluster {i+1} – {len(arcs)} selected arcs:")
    for arc in sorted(arcs):
        print(f"  {arc}")

for i, arcs in enumerate(unselected):
    print(f"\n❌ Cluster {i+1} – {len(arcs)} unresolved arcs:")
    for arc in sorted(arcs):
        print(f"  {arc}")
