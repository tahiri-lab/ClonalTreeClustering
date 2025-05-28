import matplotlib.pyplot as plt
import re

# ------------------------------
# Parse a Newick-formatted tree
# Supports named internal nodes
# ------------------------------
def parse_newick(newick):
    # Tokenize the Newick string into symbols and node names
    tokens = re.findall(r'\(|\)|,|[^(),;]+', newick)
    stack = []
    node_id = [0]

    # Helper function to create a new node with a unique ID
    def new_node(name=None):
        nid = f"n{node_id[0]}"
        node_id[0] += 1
        return {"id": nid, "name": name, "children": []}

    current = new_node("root")  # Start with a root node
    for token in tokens:
        if token == "(":
            # Start a new subtree
            child = new_node()
            current["children"].append(child)
            stack.append(current)
            current = child
        elif token == ",":
            # Start a new sibling node
            current = stack[-1]
            sibling = new_node()
            current["children"].append(sibling)
            current = sibling
        elif token == ")":
            # Close the current subtree
            current = stack.pop()
        else:
            # Assign the node name
            current["name"] = token
    return current

# ---------------------------------------------------------
# Assign (x, y) positions to nodes using post-order layout
# Leaves are positioned from left to right; parents centered
# ---------------------------------------------------------
def assign_positions(node, depth=0, x_pos=[0], positions=None):
    if positions is None:
        positions = {}

    if not node["children"]:
        # Leaf node: assign x and y based on depth and order
        x = x_pos[0]
        positions[node["id"]] = (x, -depth)
        x_pos[0] += 1
    else:
        # Internal node: assign positions to children first
        for child in node["children"]:
            assign_positions(child, depth + 1, x_pos, positions)
        # Center the parent above its children
        child_xs = [positions[child["id"]][0] for child in node["children"]]
        x = sum(child_xs) / len(child_xs)
        positions[node["id"]] = (x, -depth)
    return positions

# -------------------------------------
# Recursively draw the tree using Matplotlib
# Draws node labels and connecting lines
# -------------------------------------
def draw_tree(ax, node, positions):
    x, y = positions[node["id"]]
    # Draw the node label (blue for internal, green for leaf)
    ax.text(x, y, node["name"] or "", ha="center", va="center",
            bbox=dict(facecolor="lightblue" if node["children"] else "lightgreen", edgecolor="black"))
    # Draw edges to children
    for child in node["children"]:
        cx, cy = positions[child["id"]]
        ax.plot([x, cx], [y, cy], 'k-')  # Draw line from parent to child
        draw_tree(ax, child, positions)  # Recursive draw

# ----------------------------
# Example: draw several trees
# ----------------------------
newick_trees = [
    "(((12)2)9,((1,4)7,(5,6,13)3)10,14,(11)15,(16)8)N;",
    "(((2)9)12,((3,1)7,(13,5,6)4)10,14,(15)11,(8)16)N;",
    "(((12)2)9,((13,4)7,(5,6,1)3)10,14,(15)11,(16)8)N;",
    "(((9)2)12,((1,3)7,(4,5,13)6)10,14,(11)15,(8)16)N;",
    "(((12)9)2,((6,1)7,(13,4,5)3)10,14,(11)15,(16)8)N;"
]

# Create one subplot per tree
fig, axs = plt.subplots(len(newick_trees), 1, figsize=(12, len(newick_trees) * 4))
if len(newick_trees) == 1:
    axs = [axs]  # Ensure axs is iterable even for one tree

# Parse, position, and draw each tree
for i, (newick, ax) in enumerate(zip(newick_trees, axs)):
    tree = parse_newick(newick)
    positions = assign_positions(tree)
    draw_tree(ax, tree, positions)
    ax.set_title(f"Tree {i+1}", fontsize=14)
    ax.axis("off")  # Hide axes for clean output

plt.tight_layout()
plt.show()
