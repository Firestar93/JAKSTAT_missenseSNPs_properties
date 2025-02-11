import pandas as pd
import matplotlib.pyplot as plt
from upsetplot import UpSet, from_indicators
import matplotlib
matplotlib.use("TkAgg")  # Use a compatible backend
import warnings
warnings.simplefilter(action="ignore", category=FutureWarning)

output_file_path_table = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\STAT_JAK_HKs_paper_request\\NucleotideSequence\\EnzymeCuts\\PythonOutput\\filtered_enzyme_protein.tsv"

# Define file paths
file_paths = {
    "STAT1": "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\STAT_JAK_HKs_paper_request\\NucleotideSequence\\EnzymeCuts\\PythonOutput\\STAT1.tsv",
    "STAT3": "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\STAT_JAK_HKs_paper_request\\NucleotideSequence\\EnzymeCuts\\PythonOutput\\STAT3.tsv",
    "STAT4": "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\STAT_JAK_HKs_paper_request\\NucleotideSequence\\EnzymeCuts\\PythonOutput\\STAT4.tsv",
    "STAT5B": "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\STAT_JAK_HKs_paper_request\\NucleotideSequence\\EnzymeCuts\\PythonOutput\\STAT5B.tsv",
    "JAK2": "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\STAT_JAK_HKs_paper_request\\NucleotideSequence\\EnzymeCuts\\PythonOutput\\JAK2.tsv",
    "JAK3": "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\STAT_JAK_HKs_paper_request\\NucleotideSequence\\EnzymeCuts\\PythonOutput\\JAK3.tsv",
    "TYK2": "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\STAT_JAK_HKs_paper_request\\NucleotideSequence\\EnzymeCuts\\PythonOutput\\TYK2.tsv"
}

# Load data into a dictionary
dataframes = {name: pd.read_csv(path, sep="\t") for name, path in file_paths.items()}

# Extract enzyme occurrences across all datasets
enzyme_protein_map = {}

for protein, df in dataframes.items():
    if "Enzyme" in df.columns:
        enzymes = df["Enzyme"].dropna().unique()
        for enzyme in enzymes:
            if enzyme not in enzyme_protein_map:
                enzyme_protein_map[enzyme] = set()
            enzyme_protein_map[enzyme].add(protein)

# Convert enzyme presence across proteins into a DataFrame for UpSet plot
enzyme_protein_df = pd.DataFrame.from_dict(
    {enzyme: {protein: (protein in enzyme_protein_map[enzyme]) for protein in file_paths.keys()} for enzyme in enzyme_protein_map},
    orient="index"
)

# Ensure Boolean Type
enzyme_protein_df = enzyme_protein_df.astype(bool)

# Count occurrences of each enzyme across genes
enzyme_counts = enzyme_protein_df.sum(axis=1)  # Sum occurrences per enzyme
print("\nEnzyme occurrence counts before filtering:\n", enzyme_counts)

# Remove enzymes that appear in only **one** JAK/STAT gene
enzyme_protein_df = enzyme_protein_df[enzyme_protein_df.sum(axis=1) > 1]

# Debug: Check after filtering
print("\nEnzyme occurrence counts after filtering:\n", enzyme_protein_df.sum(axis=1))

enzyme_protein_df.to_csv(output_file_path_table, sep="\t", index=True)

# Correct usage of from_indicators
upset_data = from_indicators(enzyme_protein_df)

# Generate the UpSet plot with bigger figure & better enzyme label display
fig = plt.figure(figsize=(120, 50))  # Increase figure size to prevent overlap
upset = UpSet(
    upset_data,
    subset_size="count",
    #show_percentages=True,
    orientation="horizontal",  # Enzymes on the right, JAKs/STATs on top
    show_counts=True,          # Display count labels
    element_size=30  # Adjusts font size of labels and prevents overlap
)
upset.plot(fig=fig)

# Rotate enzyme names (X-axis) for better readability
plt.xticks(rotation=45, ha="right", fontsize=10)  # Reduce font size slightly
plt.yticks(fontsize=10)  # Reduce font size slightly

# Adjust subplot spacing to prevent clipping
plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.3)

plt.suptitle("Enzyme Restriction Site Occurrences Across JAK-STAT Proteins", fontsize=20)
plt.show(block=True)  # Keep plot open
