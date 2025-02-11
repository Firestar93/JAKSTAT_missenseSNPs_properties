from Bio.Restriction import AllEnzymes
import csv

# Define the output file
output_file = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\STAT_JAK_HKs_paper_request\\NucleotideSequence\\EnzymeCuts\\PythonOutput\\restriction_enzymes.tsv"

# Open the file and write enzyme names and recognition patterns
with open(output_file, "w", newline="") as tsvfile:
    writer = csv.writer(tsvfile, delimiter="\t")
    writer.writerow(["RestrictionEnzymeName", "RecognitionPattern"])  # Header row

    for enzyme in AllEnzymes:
        recognition_pattern = enzyme.site if enzyme.site else "Unknown"
        writer.writerow([enzyme.__name__, recognition_pattern])

print(f"TSV file '{output_file}' created successfully.")
