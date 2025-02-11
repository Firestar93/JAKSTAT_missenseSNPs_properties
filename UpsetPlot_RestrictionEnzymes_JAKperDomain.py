import pandas as pd
import matplotlib.pyplot as plt
from upsetplot import UpSet, from_indicators
import re

enzyme_domain_map_output_path = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\STAT_JAK_HKs_paper_request\\NucleotideSequence\\EnzymeCuts\\PythonOutput\\UpsetPlot_JAK_enzyme_domain_map.tsv"
upsetplot_path = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\STAT_JAK_HKs_paper_request\\NucleotideSequence\\EnzymeCuts\\PythonOutput\\UpsetPlot_JAK_enzyme_domain_map.png"

# Define file paths
jak_file_names = {
    "JAK2": "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\STAT_JAK_HKs_paper_request\\NucleotideSequence\\EnzymeCuts\\PythonOutput\\JAK2.tsv",
    "JAK3": "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\STAT_JAK_HKs_paper_request\\NucleotideSequence\\EnzymeCuts\\PythonOutput\\JAK3.tsv",
    "TYK2": "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\STAT_JAK_HKs_paper_request\\NucleotideSequence\\EnzymeCuts\\PythonOutput\\TYK2.tsv"
}

# Define STAT domain positions
# Define JAK domain positions
jak_domains = {
    "JAK2": {
        "FERM": (37, 380),
        "SH2": (386, 482),
        "Pseudokinase": (545, 809),
        "Kinase": (844, 1132)
    },
    "JAK3": {
        "FERM": (24, 356),
        "SH2": (362, 458),
        "Pseudokinase": (521, 781),
        "Kinase": (822, 1124)
    },
    "TYK2": {
        "FERM": (26, 431),
        "SH2": (452, 539),
        "Pseudokinase": (589, 875),
        "Kinase": (892, 1187)
    }
}

# Initialize dictionary to store enzyme-domain associations
enzyme_domain_map = {}

# Function to extract amino acid position from mutation names (e.g., "Val266Ile" → 266)
def extract_aa_position(mutation_name):
    match = re.search(r'\d+', mutation_name)
    return int(match.group()) if match else None


# Process each STAT file
for stat, path in jak_file_names.items():
    df = pd.read_csv(path, sep="\t")

    # Ensure 'Enzyme' column exists
    if "Enzyme" in df.columns:
        mutation_columns = [col for col in df.columns if col not in ["Enzyme", "Recognition"]]

        for _, row in df.iterrows():
            enzyme = row["Enzyme"]
            if pd.notna(enzyme):

                # Initialize enzyme entry
                if enzyme not in enzyme_domain_map:
                    enzyme_domain_map[enzyme] = {domain: False for domain in jak_domains[stat]}

                for mutation in mutation_columns:
                    base_mutation = mutation.split("_")[0]  # Remove suffixes (_FullStart, _LocalStart, etc.)
                    aa_position = extract_aa_position(base_mutation)

                    if aa_position:
                        # Identify which domain this mutation belongs to
                        for domain, (start, end) in jak_domains[stat].items():
                            if start <= aa_position <= end:
                                enzyme_domain_map[enzyme][domain] = True  # Mark presence in domain

# Convert enzyme-domain mapping into a DataFrame
enzyme_domain_df = pd.DataFrame.from_dict(enzyme_domain_map, orient="index")

# Ensure Boolean Type
enzyme_domain_df = enzyme_domain_df.astype(bool)

enzyme_domain_df.to_csv(enzyme_domain_map_output_path, sep="\t", index=True)
print(f"Enzyme-Domain map saved to: {enzyme_domain_map_output_path}")

fig = plt.figure(figsize=(40, 20), dpi=300)  # High resolution
upset = UpSet(
    from_indicators(enzyme_domain_df),
    subset_size="count",
    show_percentages=False,
    orientation="horizontal",  # Enzymes as rows, Domains as columns
    show_counts=True,
    element_size=20,
    intersection_plot_elements=10,
    sort_by="cardinality"
)
upset.plot(fig=fig)

# Adjust enzyme label font size
for ax in fig.axes:
    for label in ax.get_yticklabels():
        label.set_fontsize(8)

        # Adjust subplot spacing
plt.subplots_adjust(left=0.2, right=0.95, top=0.9, bottom=0.15)
#plt.suptitle("Enzyme Restriction Sites Across STAT Domains", fontsize=20)

plt.savefig(upsetplot_path, dpi=300, bbox_inches="tight")


print("UpSet plot saved as UpSet_JAK_Domains.png")
