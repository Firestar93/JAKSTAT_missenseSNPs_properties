import pandas as pd
import matplotlib.pyplot as plt
from upsetplot import UpSet, from_indicators
import re

enzyme_domain_map_output_path = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\STAT_JAK_HKs_paper_request\\NucleotideSequence\\EnzymeCuts\\PythonOutput\\UpsetPlot_STAT_enzyme_domain_map.tsv"
upsetplot_path = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\STAT_JAK_HKs_paper_request\\NucleotideSequence\\EnzymeCuts\\PythonOutput\\UpsetPlot_STAT_enzyme_domain_map.png"

# Define file paths
stat_file_paths = {
    "STAT1": "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\STAT_JAK_HKs_paper_request\\NucleotideSequence\\EnzymeCuts\\PythonOutput\\STAT1.tsv",
    "STAT3": "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\STAT_JAK_HKs_paper_request\\NucleotideSequence\\EnzymeCuts\\PythonOutput\\STAT3.tsv",
    "STAT4": "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\STAT_JAK_HKs_paper_request\\NucleotideSequence\\EnzymeCuts\\PythonOutput\\STAT4.tsv",
    "STAT5B": "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\STAT_JAK_HKs_paper_request\\NucleotideSequence\\EnzymeCuts\\PythonOutput\\STAT5B.tsv",
    }

# Define STAT domain positions
stat_domains = {
    "STAT1": {
        "N-terminal": (1, 121),
        "Coiled-coil": (122, 313),
        "DNA-binding": (314, 477),
        "Linker": (478, 558),
        "SH2": (559, 707),
        "TAD": (708, 750)
    },
    "STAT3": {
        "N-terminal": (1, 122),
        "Coiled-coil": (123, 318),
        "DNA-binding": (319, 484),
        "Linker": (485, 565),
        "SH2": (566, 716),
        "TAD": (717, 770)
    },
    "STAT4": {
        "N-terminal": (1, 122),
        "Coiled-coil": (123, 312),
        "DNA-binding": (313, 473),
        "Linker": (474, 554),
        "SH2": (555, 700),
        "TAD": (701, 748)
    },
    "STAT5B": {
        "N-terminal": (1, 125),
        "Coiled-coil": (126, 330),
        "DNA-binding": (331, 489),
        "Linker": (490, 574),
        "SH2": (575, 712),
        "TAD": (713, 794)
    }
}

# Initialize dictionary to store enzyme-domain associations
enzyme_domain_map = {}

# Function to extract amino acid position from mutation names (e.g., "Val266Ile" → 266)
def extract_aa_position(mutation_name):
    match = re.search(r'\d+', mutation_name)
    return int(match.group()) if match else None


# Process each STAT file
for stat, path in stat_file_paths.items():
    df = pd.read_csv(path, sep="\t")

    # Ensure 'Enzyme' column exists
    if "Enzyme" in df.columns:
        mutation_columns = [col for col in df.columns if col not in ["Enzyme", "Recognition"]]

        for _, row in df.iterrows():
            enzyme = row["Enzyme"]
            if pd.notna(enzyme):

                # Initialize enzyme entry
                if enzyme not in enzyme_domain_map:
                    enzyme_domain_map[enzyme] = {domain: False for domain in stat_domains[stat]}

                for mutation in mutation_columns:
                    base_mutation = mutation.split("_")[0]  # Remove suffixes (_FullStart, _LocalStart, etc.)
                    aa_position = extract_aa_position(base_mutation)

                    if aa_position:
                        # Identify which domain this mutation belongs to
                        for domain, (start, end) in stat_domains[stat].items():
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


print("UpSet plot saved as UpSet_STAT_Domains.png")
