#!/usr/bin/env python3

from Bio.Seq import Seq
from Bio.Restriction import AllEnzymes, RestrictionBatch
import csv

def main():
    # ------------------------------------------------------
    # 1. USER INPUTS
    # ------------------------------------------------------

    # A list of mutations, each is (mutationName, start, end).
    # Here we have two examples:
    mutations = [
        ("Arg231Trp", 691, 693),
        ("Val362Phe",  1084, 1086),
        ("Gly634Glu", 1900, 1902),
        ("Ile684Ser", 2050, 2052),
        ("Arg703Trp", 2107, 2109),
        ("Gly761Val", 2281, 2283),
        ("Ala928Val", 2782, 2784),
        ("Pro1104Ala", 3310, 3312)
    ]

    # Window size: how many bases to include upstream & downstream
    window_size = 20

    #output file dir
    output_filename = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\STAT_JAK_HKs_paper_request\\NucleotideSequence\\EnzymeCuts\\PythonOutput\\TYK2.tsv"

    # Example DNA sequence (could be read from FASTA/file).
    dna_sequence = (
        "ATGCCTCTGCGCCACTGGGGGATGGCCAGGGGCAGTAAGCCCGTTGGGGATGGAGCCCAGCCCATGGCTGCCATGGGAGGCCTGAAGGTGCTTCTGCACTGGGCTGGTCCAGGCGGCGGGGAGCCCTGGGTCACTTTCAGTGAGTCATCGCTGACAGCTGAGGAAGTCTGCATCCACATTGCACATAAAGTTGGTATCACTCCTCCTTGCTTCAATCTCTTTGCCCTCTTCGATGCTCAGGCCCAAGTCTGGTTGCCCCCAAACCACATCCTAGAGATCCCCAGAGATGCAAGCCTGATGCTATATTTCCGCATAAGGTTTTATTTCCGGAACTGGCATGGCATGAATCCTCGGGAACCGGCTGTGTACCGTTGTGGGCCCCCAGGAACCGAGGCATCCTCAGATCAGACAGCACAGGGGATGCAACTCCTGGACCCAGCCTCATTTGAGTACCTCTTTGAGCAGGGCAAGCATGAGTTTGTGAATGACGTGGCATCACTGTGGGAGCTGTCGACCGAGGAGGAGATCCACCACTTTAAGAATGAGAGCCTGGGCATGGCCTTTCTGCACCTCTGTCACCTCGCTCTCCGCCATGGCATCCCCCTGGAGGAGGTGGCCAAGAAGACCAGCTTCAAGGACTGCATCCCGCGCTCCTTCCGCCGGCATATCCGGCAGCACAGCGCCCTGACCCGGCTGCGCCTTCGGAACGTCTTCCGCAGGTTCCTGCGGGACTTCCAGCCGGGCCGACTCTCCCAGCAGATGGTCATGGTCAAATACCTAGCCACACTCGAGCGGCTGGCACCCCGCTTCGGCACAGAGCGTGTGCCCGTGTGCCACCTGAGGCTGCTGGCCCAGGCCGAGGGGGAGCCCTGCTACATCCGGGACAGTGGGGTGGCCCCTACAGACCCTGGCCCTGAGTCTGCTGCTGGGCCCCCAACCCACGAGGTGCTGGTGACAGGCACTGGTGGCATCCAGTGGTGGCCAGTAGAGGAGGAGGTGAACAAGGAGGAGGGTTCTAGTGGCAGCAGTGGCAGGAACCCCCAAGCCAGCCTGTTTGGGAAGAAGGCCAAGGCTCACAAGGCAGTCGGCCAGCCGGCAGACAGGCCGCGGGAGCCACTGTGGGCCTACTTCTGTGACTTCCGGGACATCACCCACGTGGTGCTGAAAGAGCACTGTGTCAGCATCCACCGGCAGGACAACAAGTGCCTGGAGCTGAGCTTGCCTTCCCGGGCTGCGGCGCTGTCCTTCGTGTCGCTGGTGGACGGCTATTTCCGCCTGACGGCCGACTCCAGCCACTACCTGTGCCACGAGGTGGCTCCCCCACGGCTGGTGATGAGCATCCGGGATGGGATCCACGGACCCCTGCTGGAGCCATTTGTGCAGGCCAAGCTGCGGCCCGAGGACGGCCTGTACCTCATTCACTGGAGCACCAGCCACCCCTACCGCCTGATCCTCACAGTGGCCCAGCGTAGCCAGGCACCAGACGGCATGCAGAGCTTGCGGCTCCGAAAGTTCCCCATTGAGCAGCAGGACGGGGCCTTCGTGCTGGAGGGCTGGGGCCGGTCCTTCCCCAGCGTTCGGGAACTTGGGGCTGCCTTGCAGGGCTGCTTGCTGAGGGCCGGGGATGACTGCTTCTCTCTGCGTCGCTGTTGCCTGCCCCAACCAGGAGAAACCTCCAATCTCATCATCATGCGGGGGGCTCGGGCCAGCCCCAGGACACTCAACCTCAGCCAGCTCAGCTTCCACCGGGTTGACCAGAAGGAGATCACCCAGCTGTCCCACTTGGGCCAGGGCACAAGGACCAACGTGTATGAGGGCCGCCTGCGAGTGGAGGGCAGCGGGGACCCTGAGGAGGGCAAGATGGATGACGAGGACCCCCTCGTGCCTGGCAGGGACCGTGGGCAGGAGCTACGAGTGGTGCTCAAAGTGCTGGACCCTAGTCACCATGACATCGCCCTGGCCTTCTACGAGACAGCCAGCCTCATGAGCCAGGTCTCCCACACGCACCTGGCCTTCGTGCATGGCGTCTGTGTGCGCGGCCCTGAAAATATCATGGTGACAGAGTACGTGGAGCACGGACCCCTGGATGTGTGGCTGCGGAGGGAGCGGGGCCATGTGCCCATGGCTTGGAAGATGGTGGTGGCCCAGCAGCTGGCCAGCGCCCTCAGCTACCTGGAGAACAAGAACCTGGTTCATGGTAATGTGTGTGGCCGGAACATCCTGCTGGCCCGGCTGGGGTTGGCAGAGGGCACCAGCCCCTTCATCAAGCTGAGTGATCCTGGCGTGGGCCTGGGCGCCCTCTCCAGGGAGGAGCGGGTGGAGAGGATCCCCTGGCTGGCCCCCGAATGCCTACCAGGTGGGGCCAACAGCCTAAGCACCGCCATGGACAAGTGGGGGTTTGGCGCCACCCTCCTGGAGATCTGCTTTGACGGAGAGGCCCCTCTGCAGAGCCGCAGTCCCTCCGAGAAGGAGCATTTCTACCAGAGGCAGCACCGGCTGCCCGAGCCCTCCTGCCCACAGCTGGCCACACTCACCAGCCAGTGTCTGACCTATGAGCCAACCCAGAGGCCATCATTCCGCACCATCCTGCGTGACCTCACCCGGCTGCAGCCCCACAATCTTGCTGACGTCTTGACTGTGAACCCGGACTCACCGGCGTCGGACCCTACGGTTTTCCACAAGCGCTATTTGAAAAAGATCCGAGATCTGGGCGAGGGTCACTTCGGCAAGGTCAGCTTGTACTGCTACGATCCGACCAACGACGGCACTGGCGAGATGGTGGCGGTGAAAGCCCTCAAGGCAGACTGCGGCCCCCAGCACCGCTCGGGCTGGAAGCAGGAGATTGACATTCTGCGCACGCTCTACCACGAGCACATCATCAAGTACAAGGGCTGCTGCGAGGACCAAGGCGAGAAGTCGCTGCAGCTGGTCATGGAGTACGTGCCCCTGGGCAGCCTCCGAGACTACCTGCCCCGGCACAGCATCGGGCTGGCCCAGCTGCTGCTCTTCGCCCAGCAGATCTGCGAGGGCATGGCCTATCTGCACGCGCAGCACTACATCCACCGAGACCTAGCCGCGCGCAACGTGCTGCTGGACAACGACAGGCTGGTCAAGATCGGGGACTTTGGCCTAGCCAAGGCCGTGCCCGAAGGCCACGAGTACTACCGCGTGCGCGAGGATGGGGACAGCCCCGTGTTCTGGTATGCCCCAGAGTGCCTGAAGGAGTATAAGTTCTACTATGCGTCAGATGTCTGGTCCTTCGGGGTGACCCTGTATGAGCTGCTGACGCACTGTGACTCCAGCCAGAGCCCCCCCACGAAATTCCTTGAGCTCATAGGCATTGCTCAGGGTCAGATGACAGTTCTGAGACTCACTGAGTTGCTGGAACGAGGGGAGAGGCTGCCACGGCCCGACAAATGTCCCTGTGAGGTCTATCATCTCATGAAGAACTGCTGGGAGACAGAGGCGTCCTTTCGCCCAACCTTCGAGAACCTCATACCCATTCTGAAGACAGTCCATGAGAAGTACCAAGGCCAGGCCCCTTCAGTGTTCAGCGTGTGCTGA"
    )

    # ------------------------------------------------------
    # 2. CREATE A RESTRICTION BATCH FROM ALLEnzymes
    # ------------------------------------------------------
    enzyme_batch = RestrictionBatch(AllEnzymes)

    # We'll store results in a dict keyed by enzyme_name:
    # enzyme_hits = {
    #   "EcoRI": {
    #       "recognition": "GAATTC",
    #       "Val266Ile": [ (localStart, localEnd, fullStart, fullEnd), ... ],
    #       "Example2":   [ (localStart, localEnd, fullStart, fullEnd), ... ]
    #   },
    #   "BamHI": {...},
    #    ...
    # }
    enzyme_hits = {}

    # ------------------------------------------------------
    # 3. FUNCTION TO EXTRACT A WINDOW & SEARCH
    # ------------------------------------------------------
    def analyze_mutation(dna_seq, mut_name, start_pos, end_pos, win_size):
        """
        Extract a snippet ±win_size around (start_pos..end_pos) [1-based],
        search for all restriction sites in that snippet,
        and return:
          snippet_str: the snippet sequence (string)
          results: { enzyme_name: [(loc_start, loc_end, full_start, full_end), ...], ... }
          window_start: the 1-based index where the snippet begins in the full sequence
          window_end:   the 1-based index where the snippet ends in the full sequence
        """
        seq_len = len(dna_seq)

        # 3a. Clip boundaries to avoid going out of range
        window_start = max(1, start_pos - win_size)
        window_end   = min(seq_len, end_pos + win_size)

        # 3b. Convert to Python 0-based
        py_start = window_start - 1  # inclusive
        py_end   = window_end        # exclusive

        snippet_str = dna_seq[py_start:py_end]
        # Make sure snippet is uppercase for Biopython searches
        snippet_seq = Seq(snippet_str.upper())

        # 3c. Search with the RestrictionBatch
        # analysis_dict is { enzymeObject: [cut_positions (1-based in snippet_seq)], ... }
        analysis_dict = enzyme_batch.search(snippet_seq)

        # 3d. Build results: store both local snippet coords and full-sequence coords
        results = {}
        for enzyme_obj, cut_positions in analysis_dict.items():
            if not cut_positions:
                continue  # no cut sites -> skip

            enzyme_name = enzyme_obj.__name__

            # We'll store each position as (localStart, localEnd, fullStart, fullEnd).
            # For simplicity, treat the cut as a single base (start=end).
            # If you'd rather map the entire recognition motif, you'd do more logic here.
            site_positions = []
            for cpos in cut_positions:
                local_start = cpos
                local_end   = cpos
                # Convert snippet's local position to full-sequence position:
                # snippet_seq[1] corresponds to full coord = window_start
                full_start = window_start + (cpos - 1)
                full_end   = full_start
                site_positions.append((local_start, local_end, full_start, full_end))

            results[enzyme_name] = site_positions

        return snippet_str, results, window_start, window_end

    # ------------------------------------------------------
    # 4. RUN ANALYSIS FOR EACH MUTATION & ACCUMULATE
    # ------------------------------------------------------
    for (mut_name, start, end) in mutations:
        snippet_str, hits_dict, wstart, wend = analyze_mutation(
            dna_sequence, mut_name, start, end, window_size
        )

        # hits_dict: { "EcoRI": [(locS, locE, fullS, fullE), ...], "BamHI": [...] }
        # We'll add these to `enzyme_hits` by enzyme name
        for enzyme_name, pos_list in hits_dict.items():
            # If we haven't seen this enzyme yet, initialize
            if enzyme_name not in enzyme_hits:
                # Retrieve the actual enzyme object to get the recognition site
                e_obj = [e for e in enzyme_batch if e.__name__ == enzyme_name][0]
                enzyme_hits[enzyme_name] = {
                    "recognition": str(e_obj.site) if e_obj.site else "???",
                }
            # Store the snippet-level + full-level positions under this mutation
            enzyme_hits[enzyme_name][mut_name] = pos_list

    # ------------------------------------------------------
    # 5. WRITE OUTPUT TO A TSV FILE
    # ------------------------------------------------------
    # We want columns:
    #   Enzyme, Recognition
    #   For each mutation => MutationName, LocalStarts, LocalEnds, FullStarts, FullEnds

    # 5a. Build the header row
    header_cols = ["Enzyme", "Recognition"]
    for (mut_name, _, _) in mutations:
        header_cols += [
            mut_name,
            f"{mut_name}_LocalStart",
            f"{mut_name}_LocalEnd",
            f"{mut_name}_FullStart",
            f"{mut_name}_FullEnd",
        ]



    with open(output_filename, "wt", newline="", encoding="utf-8") as outfile:
        writer = csv.writer(outfile, delimiter="\t")
        # Write the header
        writer.writerow(header_cols)

        # 5b. Sort enzymes by name for consistent output
        for enzyme_name in sorted(enzyme_hits.keys()):
            recognition_site = enzyme_hits[enzyme_name]["recognition"]
            row_data = [enzyme_name, recognition_site]

            # For each mutation in the same order, gather cut positions
            for (mut_name, _, _) in mutations:
                if mut_name in enzyme_hits[enzyme_name]:
                    # pos_list is a list of (locS, locE, fullS, fullE)
                    pos_list = enzyme_hits[enzyme_name][mut_name]

                    # Build semicolon-separated lists
                    local_starts = ";".join(str(t[0]) for t in pos_list)
                    local_ends   = ";".join(str(t[1]) for t in pos_list)
                    full_starts  = ";".join(str(t[2]) for t in pos_list)
                    full_ends    = ";".join(str(t[3]) for t in pos_list)

                    row_data += [
                        mut_name,
                        local_starts,
                        local_ends,
                        full_starts,
                        full_ends,
                    ]
                else:
                    # No cut found => put placeholder
                    row_data += [mut_name, "-", "-", "-", "-"]

            writer.writerow(row_data)

    print(f"Done! Results written to '{output_filename}'")


if __name__ == "__main__":
    main()
