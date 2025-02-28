import pandas as pd
from collections import Counter
import itertools
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import matplotlib
import numpy as np
matplotlib.use('TkAgg')  # Switch to Tkinter backend

# File paths
benign_file = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\STAT_JAK_HKs_paper_request\\AminoAcidSequence\\BENIGN_mutation_flanking_sequences.tsv"
disease_file = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\STAT_JAK_HKs_paper_request\\AminoAcidSequence\\mutation_flanking_sequences.tsv"

amino_acid_comparison_file = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\STAT_JAK_HKs_paper_request\\AminoAcidSequence\\amino_acid_comparison_file.tsv"
amino_acid_combination1_comparison= "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\STAT_JAK_HKs_paper_request\\AminoAcidSequence\\amino_acid_combination1_comparison.tsv"
amino_acid_combination2_comparison= "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\STAT_JAK_HKs_paper_request\\AminoAcidSequence\\amino_acid_combination2_comparison.tsv"
amino_acid_combination3_comparison= "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\STAT_JAK_HKs_paper_request\\AminoAcidSequence\\amino_acid_combination3_comparison.tsv"
amino_acid_combination4_comparison= "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\STAT_JAK_HKs_paper_request\\AminoAcidSequence\\amino_acid_combination4_comparison.tsv"
amino_acid_combination5_comparison= "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\STAT_JAK_HKs_paper_request\\AminoAcidSequence\\amino_acid_combination5_comparison.tsv"
amino_acid_combination6_comparison= "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\STAT_JAK_HKs_paper_request\\AminoAcidSequence\\amino_acid_combination6_comparison.tsv"

plot_dir="C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\STAT_JAK_HKs_paper_request\\AminoAcidSequence\\plots\\"

# Function to read files and extract sequences
def read_sequences(file_path):
    df = pd.read_csv(file_path, sep='\t')  # Assuming tab-separated values
    df['Upstream (AA)'] = df['Upstream (AA)'].str.split('-')
    df['Downstream (AA)'] = df['Downstream (AA)'].str.split('-')
    sequences = df['Upstream (AA)'] + df['Downstream (AA)']  # Merge upstream and downstream
    return sequences.tolist()


# Function to count amino acid occurrences
def count_amino_acids(sequences):
    aa_counter = Counter()
    for seq in sequences:
        aa_counter.update(seq)
    return aa_counter


# Function to count amino acid combinations
def count_combinations(sequences, combo_size=2):
    combo_counter = Counter()
    for seq in sequences:
        combos = list(itertools.combinations(seq, combo_size))
        combo_counter.update(combos)
    return combo_counter


# Function to plot grouped bar chart
# Function to plot grouped bar chart
def plot_grouped_bar_chart(dataframe, filename, title="", xlabel="", ylabel="", threshold=3):
    filtered_df = dataframe[(dataframe['Benign'] > threshold) | (dataframe['Disease'] > threshold)]

    fig, ax = plt.subplots(figsize=(25, 12))  # Increased height for better visibility
    indices = np.arange(len(filtered_df))
    width = 0.4

    ax.bar(indices - width / 2, filtered_df['Benign'], width=width, label='Benign', color='blue', alpha=0.6)
    ax.bar(indices + width / 2, filtered_df['Disease'], width=width, label='Disease', color='red', alpha=0.6)

    formatted_labels = ['+'.join(combo) for combo in filtered_df.index]
    ax.set_xticks(indices)
    ax.set_xticklabels(formatted_labels, rotation=90, fontsize=20, fontweight='bold')
    ax.set_xlabel(xlabel, fontsize=20, fontweight='bold')
    ax.set_ylabel(ylabel, fontsize=20, fontweight='bold')
    ax.set_title(title, fontsize=16, fontweight='bold')
    ax.legend(fontsize=30)
    ax.tick_params(axis='y', labelsize=20, width=2)  # Increase y-axis tick label size
    ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))  # Ensure only full number y-ticks

    plt.tight_layout()
    plt.savefig(filename)


# Read sequences
benign_sequences = read_sequences(benign_file)
disease_sequences = read_sequences(disease_file)

# Count amino acids
benign_aa_counts = count_amino_acids(benign_sequences)
disease_aa_counts = count_amino_acids(disease_sequences)

benign_combo1_counts = count_combinations(benign_sequences, 1)
disease_combo1_counts = count_combinations(disease_sequences, 1)

# Count amino acid combinations
benign_combo2_counts = count_combinations(benign_sequences)
disease_combo2_counts = count_combinations(disease_sequences)

benign_combo3_counts = count_combinations(benign_sequences,3)
disease_combo3_counts = count_combinations(disease_sequences,3)

benign_combo4_counts = count_combinations(benign_sequences, 4)
disease_combo4_counts = count_combinations(disease_sequences, 4)

benign_combo5_counts = count_combinations(benign_sequences,5)
disease_combo5_counts = count_combinations(disease_sequences, 5)

benign_combo6_counts = count_combinations(benign_sequences, 6)
disease_combo6_counts = count_combinations(disease_sequences, 6)

# Convert to DataFrame for comparison
aa_comparison = pd.DataFrame({'Benign': benign_aa_counts, 'Disease': disease_aa_counts}).fillna(0)
combo_comparison1 = pd.DataFrame({'Benign': benign_combo1_counts, 'Disease': disease_combo1_counts}).fillna(0)
combo_comparison2 = pd.DataFrame({'Benign': benign_combo2_counts, 'Disease': disease_combo2_counts}).fillna(0)
combo_comparison3 = pd.DataFrame({'Benign': benign_combo3_counts, 'Disease': disease_combo3_counts}).fillna(0)
combo_comparison4 = pd.DataFrame({'Benign': benign_combo4_counts, 'Disease': disease_combo4_counts}).fillna(0)
combo_comparison5 = pd.DataFrame({'Benign': benign_combo5_counts, 'Disease': disease_combo5_counts}).fillna(0)
combo_comparison6 = pd.DataFrame({'Benign': benign_combo6_counts, 'Disease': disease_combo6_counts}).fillna(0)

# Save results to TSV files
aa_comparison.to_csv(amino_acid_comparison_file, sep='\t')
combo_comparison1.to_csv(amino_acid_combination1_comparison, sep='\t')
combo_comparison2.to_csv(amino_acid_combination2_comparison, sep='\t')
combo_comparison3.to_csv(amino_acid_combination3_comparison, sep='\t')
combo_comparison4.to_csv(amino_acid_combination4_comparison, sep='\t')
combo_comparison5.to_csv(amino_acid_combination5_comparison, sep='\t')
combo_comparison6.to_csv(amino_acid_combination6_comparison, sep='\t')

plot_grouped_bar_chart(combo_comparison1,plot_dir+"combo_comparison1.png")
plot_grouped_bar_chart(combo_comparison2,plot_dir+"combo_comparison2.png")
plot_grouped_bar_chart(combo_comparison3,plot_dir+"combo_comparison3.png")
plot_grouped_bar_chart(combo_comparison4,plot_dir+"combo_comparison4.png",threshold=2)
plot_grouped_bar_chart(combo_comparison5,plot_dir+"combo_comparison5.png",threshold=1)
plot_grouped_bar_chart(combo_comparison6,plot_dir+"combo_comparison6.png",threshold=1)



