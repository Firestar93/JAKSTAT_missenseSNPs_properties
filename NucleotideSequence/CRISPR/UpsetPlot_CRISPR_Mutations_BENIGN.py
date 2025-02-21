import pandas as pd
import matplotlib.pyplot as plt
from upsetplot import UpSet, from_indicators
import matplotlib

matplotlib.use("TkAgg")  # Use a compatible backend
import warnings

warnings.simplefilter(action="ignore", category=FutureWarning)

# Save the full enzyme-mutation map
full_enzyme_mutation_map_output_path = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\STAT_JAK_HKs_paper_request\\NucleotideSequence\\CRISPR\\FoundCRISPRsites\\UpsetPlot_CRISPR_BENIGN_Variants.tsv"

high_res_output_path = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\STAT_JAK_HKs_paper_request\\NucleotideSequence\\CRISPR\\FoundCRISPRsites\\UpsetPlot_CRISPR_BENIGN_Variants.png"


# Define file paths
file_paths = {
    "STAT3": "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\STAT_JAK_HKs_paper_request\\NucleotideSequence\\CRISPR\\FoundCRISPRsites\\BENIGN\\BENIGN_STAT3.tsv",
    "STAT4": "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\STAT_JAK_HKs_paper_request\\NucleotideSequence\\CRISPR\\FoundCRISPRsites\\BENIGN\\BENIGN_STAT4.tsv",
    "STAT5B": "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\STAT_JAK_HKs_paper_request\\NucleotideSequence\\CRISPR\\FoundCRISPRsites\\BENIGN\\BENIGN_STAT5B.tsv",
    "JAK2": "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\STAT_JAK_HKs_paper_request\\NucleotideSequence\\CRISPR\\FoundCRISPRsites\\BENIGN\\BENIGN_JAK2.tsv",
    "JAK3": "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\STAT_JAK_HKs_paper_request\\NucleotideSequence\\CRISPR\\FoundCRISPRsites\\BENIGN\\BENIGN_JAK3.tsv",
    "TYK2": "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\STAT_JAK_HKs_paper_request\\NucleotideSequence\\CRISPR\\FoundCRISPRsites\\BENIGN\\BENIGN_TYK2.tsv"
}

# Initialize dictionaries
enzyme_site_counts = {}  # Counts occurrences of each enzyme restriction site
enzyme_mutation_map = {}  # Maps enzymes to associated mutation names

# Function to extract base mutation names (remove suffixes like _FullStart, _LocalStart, etc.)
def extract_base_mutation(mutation_column):
    return mutation_column.split("_")[0]  # Keep only the base mutation name

# Process each file to extract mutation-enzyme associations
for gene, path in file_paths.items():
    df = pd.read_csv(path, sep="\t")

    # Ensure 'Enzyme' column exists
    if "Cas_protein" in df.columns:
        mutation_columns = [col for col in df.columns if col not in ["Cas_protein", "Recognition"]]  # Identify mutation-related columns

        for _, row in df.iterrows():
            enzyme = row["Cas_protein"]
            if pd.notna(enzyme):  # Ensure enzyme is valid

                # Count how often each enzyme is hit
                enzyme_site_counts[enzyme] = enzyme_site_counts.get(enzyme, 0) + 1

                # Track which mutations are associated with this enzyme
                if enzyme not in enzyme_mutation_map:
                    enzyme_mutation_map[enzyme] = set()

                for mutation in mutation_columns:
                    base_mutation = extract_base_mutation(mutation)  # Extract mutation name without suffix
                    if pd.notna(row[mutation]) and row[mutation] not in ["-", "",mutation]:
                        # If mutation site has an enzyme hit (not empty or "-"), add the base mutation
                        enzyme_mutation_map[enzyme].add(base_mutation)

# Convert enzyme hit counts into a DataFrame
enzyme_counts_df = pd.DataFrame.from_dict(enzyme_site_counts, orient='index', columns=['Hit_Count'])

# Sort by most frequently hit sites
enzyme_counts_df = enzyme_counts_df.sort_values(by='Hit_Count', ascending=False)


# Convert the full enzyme-mutation map to a structured DataFrame
all_mutations = set()
for mutations in enzyme_mutation_map.values():
    all_mutations.update(mutations)

full_enzyme_mutation_df = pd.DataFrame.from_dict(
    {enzyme: {mutation: (mutation in enzyme_mutation_map[enzyme]) for mutation in all_mutations}
     for enzyme in enzyme_mutation_map},
    orient="index"
)

# Ensure Boolean Type
full_enzyme_mutation_df = full_enzyme_mutation_df.astype(bool)

# Save the full enzyme-mutation map as a TSV file
full_enzyme_mutation_df.to_csv(full_enzyme_mutation_map_output_path, sep="\t", index=True)
print(f"Full Enzyme-Mutation map saved to: {full_enzyme_mutation_map_output_path}")

# Extract top 20 most frequently hit enzyme sites
top_20_enzyme_sites = enzyme_counts_df

# Filter enzyme presence data to only include top 20 enzyme sites for visualization
top_20_enzymes = top_20_enzyme_sites.index.tolist()
top_20_enzyme_mutation_df = full_enzyme_mutation_df.loc[top_20_enzymes]

top_20_enzyme_mutation_df = top_20_enzyme_mutation_df.T

# Generate the UpSet plot with enzymes on the right and mutation names on top
fig = plt.figure(figsize=(40, 50))  # Increase figure size to prevent overlap
upset = UpSet(
    from_indicators(top_20_enzyme_mutation_df),
    subset_size="count",
    show_percentages=False,
    orientation="horizontal",  # Enzymes on the right, Mutation Names on top
    show_counts=True,  # Display count labels
    element_size=12,  # Adjusts font size of labels and prevents overlap
    intersection_plot_elements=10,
    sort_by="cardinality",
)
upset.plot(fig=fig)

for ax in fig.axes:
    for label in ax.get_yticklabels():  # Get Y labels (Enzyme Names)
        label.set_fontsize(8)  # Adjust this value (Try 8, 9, or 10)

# Rotate mutation names (X-axis) for better readability
plt.xticks(rotation=45, ha="right", fontsize=50)  # Reduce font size slightly
plt.yticks(fontsize=5)  # Reduce font size slightly

# Adjust subplot spacing to prevent clipping
plt.subplots_adjust(left=0.2, right=0.95, top=0.85, bottom=0.25)

plt.savefig(high_res_output_path, dpi=300, bbox_inches="tight")  # ✅ Ensures high DPI & no cutoff


plt.suptitle("Frequently Hit CRISPR Sites Across Mutations", fontsize=20)
plt.show(block=True)  # Keep plot open