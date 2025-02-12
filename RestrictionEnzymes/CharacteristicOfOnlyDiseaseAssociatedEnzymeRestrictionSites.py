import matplotlib.pyplot as plt
import networkx as nx
import matplotlib
matplotlib.use('TkAgg')  # Switch to Tkinter backend

# Create a graph
G = nx.Graph()

# Categories and enzymes
categories = {
    "Recognize GATC": ["DpnI", "DpnII", "Sau3AI", "Bsp143I", "Asi256I", "Lcr047I", "MalI"],
    "Methylation-sensitive": ["DpnI", "DpnII", "Sau3AI"],
    "Type IIS (Cuts outside site)": ["FaiI", "BssMI"],
    "Fixed Cutting Site": ["NdeI", "DpnII", "Sau3AI"],
    "Useful in Molecular Cloning": ["BssMI", "NdeI", "DpnI"]
}

# Add nodes and edges
for category, enzymes in categories.items():
    G.add_node(category, type="category")
    for enzyme in enzymes:
        G.add_node(enzyme, type="enzyme")
        G.add_edge(category, enzyme)

# Create node colors
color_map = []
size_map = []
for node in G:
    if G.nodes[node]["type"] == "category":
        color_map.append("lightblue")  # Categories in blue
        size_map.append(2000)  # Larger size for category nodes
    else:
        color_map.append("lightcoral")  # Enzymes in red
        size_map.append(1500)  # Smaller size for enzyme nodes

# Draw the graph
plt.figure(figsize=(12, 10))
pos = nx.spring_layout(G, seed=42, k=0.8)  # Adjust k for spacing
nx.draw(
    G, pos, with_labels=True, node_color=color_map, edge_color="gray",
    node_size=size_map, font_size=10, font_weight="bold"
)

# Title
plt.title("Common Properties of Restriction Enzymes", fontsize=14, fontweight="bold")

# Show plot
plt.show()
