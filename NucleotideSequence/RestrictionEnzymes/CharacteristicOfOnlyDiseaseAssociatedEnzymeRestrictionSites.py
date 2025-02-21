import matplotlib

matplotlib.use('TkAgg')  # Use the Tkinter backend
import matplotlib.pyplot as plt
import networkx as nx

# Create a new graph
G = nx.Graph()

# Categories and enzymes
# Feel free to adjust membership based on actual literature data!
categories = {
    # Example: These enzymes recognize a GATC-containing site
    "Recognize GATC": [
        "DpnI", "DpnII", "Sau3AI", "Bsp143I", "Asi256I", "Lcr047I", "MaII"
    ],

    # Example: Methylation-sensitive or dependent
    "Methylation-sensitive": [
        "DpnI", "DpnII", "Sau3AI", "MspJI"
    ],

    # Example: Type IIS (enzymes that cut outside their recognition site)
    "Type IIS (Cuts outside site)": [
        "FaiI", "BssMI"
    ],

    # Example: These enzymes cut at a fixed site (not outside)
    "Fixed Cutting Site": [
        "NdeI", "DpnII", "Sau3AI"
    ],

    # Example: Commonly used in molecular cloning
    "Useful in Molecular Cloning": [
        "BssMI", "NdeI", "DpnI"
    ]
}

# ---------------------------------------------------------------------
# 1) Add nodes for categories (type="category") and enzymes (type="enzyme").
# 2) Add edges from each category to every enzyme that belongs to it.
# ---------------------------------------------------------------------
for category, enzymes in categories.items():
    # Add the category as a node
    G.add_node(category, type="category")

    # Add each enzyme as a node, then link it to the category
    for enzyme in enzymes:
        G.add_node(enzyme, type="enzyme")
        G.add_edge(category, enzyme)

# ---------------------------------------------------------------------
# Create a color and size map so categories and enzymes appear differently.
# ---------------------------------------------------------------------
node_colors = []
node_sizes = []
for node in G.nodes():
    if G.nodes[node]["type"] == "category":
        node_colors.append("lightblue")  # Categories in blue
        node_sizes.append(2000)  # Larger size for category nodes
    else:
        node_colors.append("lightcoral")  # Enzymes in red
        node_sizes.append(1500)  # Smaller size for enzyme nodes

# ---------------------------------------------------------------------
# Draw the graph
# ---------------------------------------------------------------------
plt.figure(figsize=(12, 10))

# Spring layout for a clean look. Adjust 'k' or 'iterations' if needed.
pos = nx.spring_layout(G, seed=42, k=0.8)

nx.draw(
    G,
    pos,
    with_labels=True,
    labels={n: n for n in G.nodes()},
    node_color=node_colors,
    node_size=node_sizes,
    edge_color="gray",
    font_size=10,
    font_weight="bold"
)

# Title
plt.title("Common Properties of Restriction Enzymes", fontsize=14, fontweight="bold")

# Display the plot
plt.show()
