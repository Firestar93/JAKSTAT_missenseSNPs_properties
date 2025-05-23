import csv
from collections import Counter
import re
from Bio.Data import CodonTable

# Output file path
output_file = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\STAT_JAK_HKs_paper_request\\AminoAcidSequence\\BENIGN_mutation_flanking_sequences.tsv"

gene_data = {
    "STAT3": {
        "sequence": (
            "ATGGCCCAATGGAATCAGCTACAGCAGCTTGACACACGGTACCTGGAGCAGCTCCATCAGCTCTACAGTGACAGCTTCCCAATGGAGCTGCGGCAGTTTCTGGCCCCTTGGATTGAGAGTCAAGATTGGGCATATGCGGCCAGCAAAGAATCACATGCCACTTTGGTGTTTCATAATCTCCTGGGAGAGATTGACCAGCAGTATAGCCGCTTCCTGCAAGAGTCGAATGTTCTCTATCAGCACAATCTACGAAGAATCAAGCAGTTTCTTCAGAGCAGGTATCTTGAGAAGCCAATGGAGATTGCCCGGATTGTGGCCCGGTGCCTGTGGGAAGAATCACGCCTTCTACAGACTGCAGCCACTGCGGCCCAGCAAGGGGGCCAGGCCAACCACCCCACAGCAGCCGTGGTGACGGAGAAGCAGCAGATGCTGGAGCAGCACCTTCAGGATGTCCGGAAGAGAGTGCAGGATCTAGAACAGAAAATGAAAGTGGTAGAGAATCTCCAGGATGACTTTGATTTCAACTATAAAACCCTCAAGAGTCAAGGAGACATGCAAGATCTGAATGGAAACAACCAGTCAGTGACCAGGCAGAAGATGCAGCAGCTGGAACAGATGCTCACTGCGCTGGACCAGATGCGGAGAAGCATCGTGAGTGAGCTGGCGGGGCTTTTGTCAGCGATGGAGTACGTGCAGAAAACTCTCACGGACGAGGAGCTGGCTGACTGGAAGAGGCGGCAACAGATTGCCTGCATTGGAGGCCCGCCCAACATCTGCCTAGATCGGCTAGAAAACTGGATAACGTCATTAGCAGAATCTCAACTTCAGACCCGTCAACAAATTAAGAAACTGGAGGAGTTGCAGCAAAAAGTTTCCTACAAAGGGGACCCCATTGTACAGCACCGGCCGATGCTGGAGGAGAGAATCGTGGAGCTGTTTAGAAACTTAATGAAAAGTGCCTTTGTGGTGGAGCGGCAGCCCTGCATGCCCATGCATCCTGACCGGCCCCTCGTCATCAAGACCGGCGTCCAGTTCACTACTAAAGTCAGGTTGCTGGTCAAATTCCCTGAGTTGAATTATCAGCTTAAAATTAAAGTGTGCATTGACAAAGACTCTGGGGACGTTGCAGCTCTCAGAGGATCCCGGAAATTTAACATTCTGGGCACAAACACAAAAGTGATGAACATGGAAGAATCCAACAACGGCAGCCTCTCTGCAGAATTCAAACACTTGACCCTGAGGGAGCAGAGATGTGGGAATGGGGGCCGAGCCAATTGTGATGCTTCCCTGATTGTGACTGAGGAGCTGCACCTGATCACCTTTGAGACCGAGGTGTATCACCAAGGCCTCAAGATTGACCTAGAGACCCACTCCTTGCCAGTTGTGGTGATCTCCAACATCTGTCAGATGCCAAATGCCTGGGCGTCCATCCTGTGGTACAACATGCTGACCAACAATCCCAAGAATGTAAACTTTTTTACCAAGCCCCCAATTGGAACCTGGGATCAAGTGGCCGAGGTCCTGAGCTGGCAGTTCTCCTCCACCACCAAGCGAGGACTGAGCATCGAGCAGCTGACTACACTGGCAGAGAAACTCTTGGGACCTGGTGTGAATTATTCAGGGTGTCAGATCACATGGGCTAAATTTTGCAAAGAAAACATGGCTGGCAAGGGCTTCTCCTTCTGGGTCTGGCTGGACAATATCATTGACCTTGTGAAAAAGTACATCCTGGCCCTTTGGAACGAAGGGTACATCATGGGCTTTATCAGTAAGGAGCGGGAGCGGGCCATCTTGAGCACTAAGCCTCCAGGCACCTTCCTGCTAAGATTCAGTGAAAGCAGCAAAGAAGGAGGCGTCACTTTCACTTGGGTGGAGAAGGACATCAGCGGTAAGACCCAGATCCAGTCCGTGGAACCATACACAAAGCAGCAGCTGAACAACATGTCATTTGCTGAAATCATCATGGGCTATAAGATCATGGATGCTACCAATATCCTGGTGTCTCCACTGGTCTATCTCTATCCTGACATTCCCAAGGAGGAGGCATTCGGAAAGTATTGTCGGCCAGAGAGCCAGGAGCATCCTGAAGCTGACCCAGGTAGCGCTGCCCCATACCTGAAGACCAAGTTTATCTGTGTGACACCAACGACCTGCAGCAATACCATTGACCTGCCGATGTCCCCCCGCACTTTAGATTCATTGATGCAGTTTGGAAATAATGGTGAAGGTGCTGAACCCTCAGCAGGAGGGCAGTTTGAGTCCCTCACCTTTGACATGGAGTTGACCTCGGAGTGCGCTACCTCCCCCATGTGA"
        # SHORT placeholder
        ),
        "mutations": [
            ("Val461Leu", 1381, 1383), ("Ile498Val", 1492,1494), ("Val507Phe", 1519, 1521), ("Ala702Thr", 2104, 2106), ("Ser763Leu", 2287, 2289),
        ],
    },
    "STAT4": {
        "sequence": (
            "ATGTCTCAGTGGAATCAAGTCCAACAGTTAGAAATCAAGTTTTTGGAGCAGGTGGATCAATTCTATGATGACAACTTTCCCATGGAAATTCGGCATCTGTTGGCCCAATGGATTGAAAATCAAGACTGGGAGGCAGCTTCTAACAATGAAACCATGGCAACGATTCTTCTTCAAAACTTGTTAATACAACTGGATGAACAGTTAGGTCGTGTTTCCAAAGAGAAAAACCTACTCTTGATACACAATCTAAAAAGAATTAGGAAGGTCCTTCAGGGAAAATTTCATGGAAATCCAATGCATGTAGCTGTGGTTATTTCAAACTGTTTAAGGGAAGAGAGGAGAATATTGGCTGCAGCCAACATGCCTGTCCAGGGGCCTCTAGAGAAATCCTTACAAAGTTCTTCAGTTTCAGAAAGACAGAGGAATGTGGAGCACAAAGTGGCTGCCATTAAAAACAGTGTGCAGATGACAGAACAAGATACCAAATACTTAGAAGATCTGCAAGACGAATTTGACTACAGGTATAAAACAATTCAGACAATGGATCAGAGTGACAAGAATAGTGCCATGGTGAATCAGGAAGTTTTGACACTGCAGGAAATGCTTAACAGCCTCGATTTCAAGAGAAAGGAGGCTCTCAGTAAAATGACCCAAATCATCCATGAGACAGACCTGTTAATGAACACCATGCTCATAGAAGAGCTGCAAGACTGGAAGCGGCGGCAGCAAATCGCCTGCATCGGGGGTCCACTCCACAATGGGCTCGACCAGCTTCAGAACTGCTTTACACTATTGGCAGAAAGTCTTTTCCAACTGAGAAGGCAATTGGAGAAACTAGAGGAGCAATCTACCAAAATGACATATGAAGGTGATCCCATTCCAATGCAAAGAACTCACATGCTAGAAAGAGTCACCTTCTTGATCTACAACCTTTTCAAGAACTCATTTGTGGTTGAGCGACAGCCATGTATGCCAACCCACCCTCAGAGGCCGTTGGTACTTAAAACCCTAATTCAGTTCACTGTAAAACTAAGGCTACTAATAAAATTGCCAGAACTAAACTATCAGGTAAAGGTTAAGGCATCAATTGACAAGAATGTTTCAACTCTAAGCAACCGAAGATTTGTACTTTGTGGAACTAATGTCAAAGCCATGTCTATTGAAGAATCTTCCAATGGGAGTCTCTCAGTAGAATTTCGACATTTGCAACCAAAGGAAATGAAGTCCAGTGCTGGAGGTAAAGGAAATGAGGGCTGTCACATGGTGACTGAAGAACTTCATTCCATAACGTTTGAAACACAGATCTGCCTCTATGGCCTGACCATAGATTTGGAGACCAGCTCATTGCCTGTGGTGATGATTTCCAATGTCAGTCAGTTACCTAATGCTTGGGCATCCATCATTTGGTACAACGTGTCAACCAACGATTCCCAGAACTTGGTTTTCTTTAATAATCCTCCACCTGCCACATTGAGTCAACTACTGGAGGTGATGAGCTGGCAGTTTTCATCGTACGTTGGTCGTGGTCTTAACTCAGATCAACTCCATATGCTGGCAGAGAAGCTTACAGTCCAATCTAGCTACAGTGATGGTCACCTCACCTGGGCCAAGTTCTGCAAGGAACATTTACCTGGTAAATCATTTACCTTTTGGACATGGCTTGAAGCAATATTGGATCTAATTAAGAAACACATTCTTCCCCTTTGGATTGATGGGTATGTCATGGGCTTTGTTAGCAAAGAGAAGGAACGGCTGTTGCTAAAGGATAAAATGCCTGGCACCTTTTTATTAAGATTCAGTGAAAGCCATCTCGGAGGAATAACTTTCACCTGGGTGGACCATTCTGAAAGTGGGGAAGTGAGATTCCACTCTGTAGAACCCTACAATAAAGGCCGGTTGTCTGCTCTGCCATTCGCTGACATCCTGCGAGACTACAAAGTTATTATGGCTGAAAACATTCCTGAAAACCCTCTGAAGTACCTATATCCTGACATTCCCAAAGACAAAGCCTTCGGTAAACACTACAGCTCTCAGCCTTGCGAAGTTTCAAGACCAACAGAAAGGGGTGACAAAGGTTATGTTCCTTCTGTTTTTATCCCCATCTCAACAATCCGAAGTGATTCAACAGAGCCACATTCTCCATCAGACCTTCTTCCCATGTCTCCAAGTGTGTATGCGGTGTTGAGAGAAAACCTGAGTCCCACAACAATTGAAACTGCAATGAAGTCTCCTTATTCTGCTGAATGA"
        ),
        "mutations": [
            ("Ile115Val", 343, 345), ("Arg240Gln", 718, 720), ("Leu269Ile", 805, 807), ("Thr298Ile", 892, 894),
        ],
    },
    "STAT5B": {
        "sequence": (
            "ATGGCTGTGTGGATACAAGCTCAGCAGCTCCAAGGAGAAGCCCTTCATCAGATGCAAGCGTTATATGGCCAGCATTTTCCCATTGAGGTGCGGCATTATTTATCCCAGTGGATTGAAAGCCAAGCATGGGACTCAGTAGATCTTGATAATCCACAGGAGAACATTAAGGCCACCCAGCTCCTGGAGGGCCTGGTGCAGGAGCTGCAGAAGAAGGCAGAGCACCAGGTGGGGGAAGATGGGTTTTTACTGAAGATCAAGCTGGGGCACTATGCCACACAGCTCCAGAACACGTATGACCGCTGCCCCATGGAGCTGGTCCGCTGCATCCGCCATATATTGTACAATGAACAGAGGTTGGTCCGAGAAGCCAACAATGGTAGCTCTCCAGCTGGAAGCCTTGCTGATGCCATGTCCCAGAAACACCTCCAGATCAACCAGACGTTTGAGGAGCTGCGACTGGTCACGCAGGACACAGAGAATGAGTTAAAAAAGCTGCAGCAGACTCAGGAGTACTTCATCATCCAGTACCAGGAGAGCCTGAGGATCCAAGCTCAGTTTGGCCCGCTGGCCCAGCTGAGCCCCCAGGAGCGTCTGAGCCGGGAGACGGCCCTCCAGCAGAAGCAGGTGTCTCTGGAGGCCTGGTTGCAGCGTGAGGCACAGACACTGCAGCAGTACCGCGTGGAGCTGGCCGAGAAGCACCAGAAGACCCTGCAGCTGCTGCGGAAGCAGCAGACCATCATCCTGGATGACGAGCTGATCCAGTGGAAGCGGCGGCAGCAGCTGGCCGGGAACGGCGGGCCCCCCGAGGGCAGCCTGGACGTGCTACAGTCCTGGTGTGAGAAGTTGGCCGAGATCATCTGGCAGAACCGGCAGCAGATCCGCAGGGCTGAGCACCTCTGCCAGCAGCTGCCCATCCCCGGCCCAGTGGAGGAGATGCTGGCCGAGGTCAACGCCACCATCACGGACATTATCTCAGCCCTGGTGACCAGCACGTTCATCATTGAGAAGCAGCCTCCTCAGGTCCTGAAGACCCAGACCAAGTTTGCAGCCACTGTGCGCCTGCTGGTGGGCGGGAAGCTGAACGTGCACATGAACCCCCCCCAGGTGAAGGCCACCATCATCAGTGAGCAGCAGGCCAAGTCTCTGCTCAAGAACGAGAACACCCGCAATGATTACAGTGGCGAGATCTTGAACAACTGCTGCGTCATGGAGTACCACCAAGCCACAGGCACCCTTAGTGCCCACTTCAGGAATATGTCCCTGAAACGAATTAAGAGGTCAGACCGTCGTGGGGCAGAGTCGGTGACAGAAGAAAAATTTACAATCCTGTTTGAATCCCAGTTCAGTGTTGGTGGAAATGAGCTGGTTTTTCAAGTCAAGACCCTGTCCCTGCCAGTGGTGGTGATCGTTCATGGCAGCCAGGACAACAATGCGACGGCCACTGTTCTCTGGGACAATGCTTTTGCAGAGCCTGGCAGGGTGCCATTTGCCGTGCCTGACAAAGTGCTGTGGCCACAGCTGTGTGAGGCGCTCAACATGAAATTCAAGGCCGAAGTGCAGAGCAACCGGGGCCTGACCAAGGAGAACCTCGTGTTCCTGGCGCAGAAACTGTTCAACAACAGCAGCAGCCACCTGGAGGACTACAGTGGCCTGTCTGTGTCCTGGTCCCAGTTCAACAGGGAGAATTTACCAGGACGGAATTACACTTTCTGGCAATGGTTTGACGGTGTGATGGAAGTGTTAAAAAAACATCTCAAGCCTCATTGGAATGATGGGGCCATTTTGGGGTTTGTAAACAAGCAACAGGCCCATGACCTACTCATTAACAAGCCAGATGGGACCTTCCTCCTGAGATTCAGTGACTCAGAAATTGGCGGCATCACCATTGCTTGGAAGTTTGATTCTCAGGAAAGAATGTTTTGGAATCTGATGCCTTTTACCACCAGAGACTTCTCCATTCGGTCCCTAGCCGACCGCTTGGGAGACTTGAATTACCTTATCTACGTGTTTCCTGATCGGCCAAAAGATGAAGTATACTCCAAATACTACACACCAGTTCCCTGCGAGTCTGCTACTGCTAAAGCTGTTGATGGATACGTGAAGCCACAGATCAAGCAAGTGGTCCCTGAGTTTGTGAACGCATCTGCAGATGCCGGGGGCGGCAGCGCCACGTACATGGACCAGGCCCCCTCCCCAGCTGTGTGTCCCCAGGCTCACTATAACATGTACCCACAGAACCCTGACTCAGTCCTTGACACCGATGGGGACTTCGATCTGGAGGACACAATGGACGTAGCGCGGCGTGTGGAGGAGCTCCTGGGCCGGCCAATGGACAGTCAGTGGATCCCGCACGCACAATCGTGA"
        ),
        "mutations": [
            ("Ala130Val", 388, 390), ("Glu315Ala", 943, 945),
        ],
    },
    "JAK2": {
        "sequence": (
            "ATGGGAATGGCCTGCCTTACGATGACAGAAATGGAGGGAACATCCACCTCTTCTATATATCAGAATGGTGATATTTCTGGAAATGCCAATTCTATGAAGCAAATAGATCCAGTTCTTCAGGTGTATCTTTACCATTCCCTTGGGAAATCTGAGGCAGATTATCTGACCTTTCCATCTGGGGAGTATGTTGCAGAAGAAATCTGTATTGCTGCTTCTAAAGCTTGTGGTATCACACCTGTGTATCATAATATGTTTGCTTTAATGAGTGAAACAGAAAGGATCTGGTATCCACCCAACCATGTCTTCCATATAGATGAGTCAACCAGGCATAATGTACTCTACAGAATAAGATTTTACTTTCCTCGTTGGTATTGCAGTGGCAGCAACAGAGCCTATCGGCATGGAATATCTCGAGGTGCTGAAGCTCCTCTTCTTGATGACTTTGTCATGTCTTACCTCTTTGCTCAGTGGCGGCATGATTTTGTGCACGGATGGATAAAAGTACCTGTGACTCATGAAACACAGGAAGAATGTCTTGGGATGGCAGTGTTAGATATGATGAGAATAGCCAAAGAAAACGATCAAACCCCACTGGCCATCTATAACTCTATCAGCTACAAGACATTCTTACCAAAATGTATTCGAGCAAAGATCCAAGACTATCATATTTTGACAAGGAAGCGAATAAGGTACAGATTTCGCAGATTTATTCAGCAATTCAGCCAATGCAAAGCCACTGCCAGAAACTTGAAACTTAAGTATCTTATAAATCTGGAAACTCTGCAGTCTGCCTTCTACACAGAGAAATTTGAAGTAAAAGAACCTGGAAGTGGTCCTTCAGGTGAGGAGATTTTTGCAACCATTATAATAACTGGAAACGGTGGAATTCAGTGGTCAAGAGGGAAACATAAAGAAAGTGAGACACTGACAGAACAGGATTTACAGTTATATTGCGATTTTCCTAATATTATTGATGTCAGTATTAAGCAAGCAAACCAAGAGGGTTCAAATGAAAGCCGAGTTGTAACTATCCATAAGCAAGATGGTAAAAATCTGGAAATTGAACTTAGCTCATTAAGGGAAGCTTTGTCTTTCGTGTCATTAATTGATGGATATTATAGATTAACTGCAGATGCACATCATTACCTCTGTAAAGAAGTAGCACCTCCAGCCGTGCTTGAAAATATACAAAGCAACTGTCATGGCCCAATTTCGATGGATTTTGCCATTAGTAAACTGAAGAAAGCAGGTAATCAGACTGGACTGTATGTACTTCGATGCAGTCCTAAGGACTTTAATAAATATTTTTTGACTTTTGCTGTCGAGCGAGAAAATGTCATTGAATATAAACACTGTTTGATTACAAAAAATGAGAATGAAGAGTACAACCTCAGTGGGACAAAGAAGAACTTCAGCAGTCTTAAAGATCTTTTGAATTGTTACCAGATGGAAACTGTTCGCTCAGACAATATAATTTTCCAGTTTACTAAATGCTGTCCCCCAAAGCCAAAAGATAAATCAAACCTTCTAGTCTTCAGAACGAATGGTGTTTCTGATGTACCAACCTCACCAACATTACAGAGGCCTACTCATATGAACCAAATGGTGTTTCACAAAATCAGAAATGAAGATTTGATATTTAATGAAAGCCTTGGCCAAGGCACTTTTACAAAGATTTTTAAAGGCGTACGAAGAGAAGTAGGAGACTACGGTCAACTGCATGAAACAGAAGTTCTTTTAAAAGTTCTGGATAAAGCACACAGAAACTATTCAGAGTCTTTCTTTGAAGCAGCAAGTATGATGAGCAAGCTTTCTCACAAGCATTTGGTTTTAAATTATGGAGTATGTGTCTGTGGAGACGAGAATATTCTGGTTCAGGAGTTTGTAAAATTTGGATCACTAGATACATATCTGAAAAAGAATAAAAATTGTATAAATATATTATGGAAACTTGAAGTTGCTAAACAGTTGGCATGGGCCATGCATTTTCTAGAAGAAAACACCCTTATTCATGGGAATGTATGTGCCAAAAATATTCTGCTTATCAGAGAAGAAGACAGGAAGACAGGAAATCCTCCTTTCATCAAACTTAGTGATCCTGGCATTAGTATTACAGTTTTGCCAAAGGACATTCTTCAGGAGAGAATACCATGGGTACCACCTGAATGCATTGAAAATCCTAAAAATTTAAATTTGGCAACAGACAAATGGAGTTTTGGTACCACTTTGTGGGAAATCTGCAGTGGAGGAGATAAACCTCTAAGTGCTCTGGATTCTCAAAGAAAGCTACAATTTTATGAAGATAGGCATCAGCTTCCTGCACCAAAGTGGGCAGAATTAGCAAACCTTATAAATAATTGTATGGATTATGAACCAGATTTCAGGCCTTCTTTCAGAGCCATCATACGAGATCTTAACAGTTTGTTTACTCCAGATTATGAACTATTAACAGAAAATGACATGTTACCAAATATGAGGATAGGTGCCCTGGGGTTTTCTGGTGCCTTTGAAGACCGGGATCCTACACAGTTTGAAGAGAGACATTTGAAATTTCTACAGCAACTTGGCAAGGGTAATTTTGGGAGTGTGGAGATGTGCCGGTATGACCCTCTACAGGACAACACTGGGGAGGTGGTCGCTGTAAAAAAGCTTCAGCATAGTACTGAAGAGCACCTAAGAGACTTTGAAAGGGAAATTGAAATCCTGAAATCCCTACAGCATGACAACATTGTAAAGTACAAGGGAGTGTGCTACAGTGCTGGTCGGCGTAATCTAAAATTAATTATGGAATATTTACCATATGGAAGTTTACGAGACTATCTTCAAAAACATAAAGAACGGATAGATCACATAAAACTTCTGCAGTACACATCTCAGATATGCAAGGGTATGGAGTATCTTGGTACAAAAAGGTATATCCACAGGGATCTGGCAACGAGAAATATATTGGTGGAGAACGAGAACAGAGTTAAAATTGGAGATTTTGGGTTAACCAAAGTCTTGCCACAAGACAAAGAATACTATAAAGTAAAAGAACCTGGTGAAAGTCCCATATTCTGGTATGCTCCAGAATCACTGACAGAGAGCAAGTTTTCTGTGGCCTCAGATGTTTGGAGCTTTGGAGTGGTTCTGTATGAACTTTTCACATACATTGAGAAGAGTAAAAGTCCACCAGCGGAATTTATGCGTATGATTGGCAATGACAAACAAGGACAGATGATCGTGTTCCATTTGATAGAACTTTTGAAGAATAATGGAAGATTACCAAGACCAGATGGATGCCCAGATGAGATCTATATGATCATGACAGAATGCTGGAACAATAATGTAAATCAACGCCCCTCCTTTAGGGATCTAGCTCTTCGAGTGGATCAAATAAGGGATAACATGGCTGGATGA"
        ),
        "mutations": [
            ("Leu113Val", 337 , 339), ("Gly127Asp", 379, 381), ("Lys244Arg", 730, 732), ("Asn337Asp", 1009, 1011), ("Val392Met", 1174, 1176),
            ("Cys480Phe", 1438, 1440), ("Leu892Val", 2674, 2676),
        ],
    },
    "JAK3": {
        "sequence": (
            "ATGGCACCTCCAAGTGAAGAGACGCCCCTGATCCCTCAGCGTTCATGCAGCCTCTTGTCCACGGAGGCTGGTGCCCTGCATGTGCTGCTGCCCGCTCGGGGCCCCGGGCCCCCCCAGCGCCTATCTTTCTCCTTTGGGGACCACTTGGCTGAGGACCTGTGCGTGCAGGCTGCCAAGGCCAGCGGCATCCTGCCTGTGTACCACTCCCTCTTTGCTCTGGCCACGGAGGACCTGTCCTGCTGGTTCCCCCCGAGCCACATCTTCTCCGTGGAGGATGCCAGCACCCAAGTCCTGCTGTACAGGATTCGCTTTTACTTCCCCAATTGGTTTGGGCTGGAGAAGTGCCACCGCTTCGGGCTACGCAAGGATTTGGCCAGTGCTATCCTTGACCTGCCAGTCCTGGAGCACCTCTTTGCCCAGCACCGCAGTGACCTGGTGAGTGGGCGCCTCCCCGTGGGCCTCAGTCTCAAGGAGCAGGGTGAGTGTCTCAGCCTGGCCGTGTTGGACCTGGCCCGGATGGCGCGAGAGCAGGCCCAGCGGCCGGGAGAGCTGCTGAAGACTGTCAGCTACAAGGCCTGCCTACCCCCAAGCCTGCGCGACCTGATCCAGGGCCTGAGCTTCGTGACGCGGAGGCGTATTCGGAGGACGGTGCGCAGAGCCCTGCGCCGCGTGGCCGCCTGCCAGGCAGACCGGCACTCGCTCATGGCCAAGTACATCATGGACCTGGAGCGGCTGGATCCAGCCGGGGCCGCCGAGACCTTCCACGTGGGCCTCCCTGGGGCCCTTGGTGGCCACGACGGGCTGGGGCTGCTCCGCGTGGCTGGTGACGGCGGCATCGCCTGGACCCAGGGAGAACAGGAGGTCCTCCAGCCCTTCTGCGACTTTCCAGAAATCGTAGACATTAGCATCAAGCAGGCCCCGCGCGTTGGCCCGGCCGGAGAGCACCGCCTGGTCACTGTTACCAGGACAGACAACCAGATTTTAGAGGCCGAGTTCCCAGGGCTGCCCGAGGCTCTGTCGTTCGTGGCGCTCGTGGACGGCTACTTCCGGCTGACCACGGACTCCCAGCACTTCTTCTGCAAGGAGGTGGCACCGCCGAGGCTGCTGGAGGAAGTGGCCGAGCAGTGCCACGGCCCCATCACTCTGGACTTTGCCATCAACAAGCTCAAGACTGGGGGCTCACGTCCTGGCTCCTATGTTCTCCGCCGCAGCCCCCAGGACTTTGACAGCTTCCTCCTCACTGTCTGTGTCCAGAACCCCCTTGGTCCTGATTATAAGGGCTGCCTCATCCGGCGCAGCCCCACAGGAACCTTCCTTCTGGTTGGCCTCAGCCGACCCCACAGCAGTCTTCGAGAGCTCCTGGCAACCTGCTGGGATGGGGGGCTGCACGTAGATGGGGTGGCAGTGACCCTCACTTCCTGCTGTATCCCCAGACCCAAAGAAAAGTCCAACCTGATCGTGGTCCAGAGAGGTCACAGCCCACCCACATCATCCTTGGTTCAGCCCCAATCCCAATACCAGCTGAGTCAGATGACATTTCACAAGATCCCTGCTGACAGCCTGGAGTGGCATGAGAACCTGGGCCATGGGTCCTTCACCAAGATTTACCGGGGCTGTCGCCATGAGGTGGTGGATGGGGAGGCCCGAAAGACAGAGGTGCTGCTGAAGGTCATGGATGCCAAGCACAAGAACTGCATGGAGTCATTCCTGGAAGCAGCGAGCTTGATGAGCCAAGTGTCGTACCGGCATCTCGTGCTGCTCCACGGCGTGTGCATGGCTGGAGACAGCACCATGGTGCAGGAATTTGTACACCTGGGGGCCATAGACATGTATCTGCGAAAACGTGGCCACCTGGTGCCAGCCAGCTGGAAGCTGCAGGTGGTCAAACAGCTGGCCTACGCCCTCAACTATCTGGAGGACAAAGGCCTGCCCCATGGCAATGTCTCTGCCCGGAAGGTGCTCCTGGCTCGGGAGGGGGCTGATGGGAGCCCGCCCTTCATCAAGCTGAGTGACCCTGGGGTCAGCCCCGCTGTGTTAAGCCTGGAGATGCTCACCGACAGGATCCCCTGGGTGGCCCCCGAGTGTCTCCGGGAGGCGCAGACACTTAGCTTGGAAGCTGACAAGTGGGGCTTCGGCGCCACGGTCTGGGAAGTGTTTAGTGGCGTCACCATGCCCATCAGTGCCCTGGATCCTGCTAAGAAACTCCAATTTTATGAGGACCGGCAGCAGCTGCCGGCCCCCAAGTGGACAGAGCTGGCCCTGCTGATTCAACAGTGCATGGCCTATGAGCCGGTCCAGAGGCCCTCCTTCCGAGCCGTCATTCGTGACCTCAATAGCCTCATCTCTTCAGACTATGAGCTCCTCTCAGACCCCACACCTGGTGCCCTGGCACCTCGTGATGGGCTGTGGAATGGTGCCCAGCTCTATGCCTGCCAAGACCCCACGATCTTCGAGGAGAGACACCTCAAGTACATCTCACAGCTGGGCAAGGGCAACTTTGGCAGCGTGGAGCTGTGCCGCTATGACCCGCTAGGCGACAATACAGGTGCCCTGGTGGCCGTGAAACAGCTGCAGCACAGCGGGCCAGACCAGCAGAGGGACTTTCAGCGGGAGATTCAGATCCTCAAAGCACTGCACAGTGATTTCATTGTCAAGTATCGTGGTGTCAGCTATGGCCCGGGCCGCCAGAGCCTGCGGCTGGTCATGGAGTACCTGCCCAGCGGCTGCTTGCGCGACTTCCTGCAGCGGCACCGCGCGCGCCTCGATGCCAGCCGCCTCCTTCTCTATTCCTCGCAGATCTGCAAGGGCATGGAGTACCTGGGCTCCCGCCGCTGCGTGCACCGCGACCTGGCCGCCCGAAACATCCTCGTGGAGAGCGAGGCACACGTCAAGATCGCTGACTTCGGCCTAGCTAAGCTGCTGCCGCTTGACAAAGACTACTACGTGGTCCGCGAGCCAGGCCAGAGCCCCATTTTCTGGTATGCCCCCGAATCCCTCTCGGACAACATCTTCTCTCGCCAGTCAGACGTCTGGAGCTTCGGGGTCGTCCTGTACGAGCTCTTCACCTACTGCGACAAAAGCTGCAGCCCCTCGGCCGAGTTCCTGCGGATGATGGGATGTGAGCGGGATGTCCCCGCCCTCTGCCGCCTCTTGGAACTGCTGGAGGAGGGCCAGAGGCTGCCGGCGCCTCCTGCCTGCCCTGCTGAGGTTCACGAGCTCATGAAGCTGTGCTGGGCCCCTAGCCCACAGGACCGGCCATCATTCAGCGCCCTGGGCCCCCAGCTGGACATGCTGTGGAGCGGAAGCCGGGGGTGTGAGACTCATGCCTTCACTGCTCACCCAGAGGGCAAACACCACTCCCTGTCCTTTTCATAG"
        ),
        "mutations": [
            ("Ile63Val", 187, 189), ("Arg121His", 361, 363), ("Val217Met", 649, 651),
            ("Pro396Leu", 1186, 1188), ("Arg451Gln", 1351, 1353), ("His879Arg", 2635, 2637), ("Ile1003Val", 3007, 3009),
            ("Ala1090Thr", 3268, 3270),
        ],
    },
    "TYK2": {
        "sequence": (
            "ATGCCTCTGCGCCACTGGGGGATGGCCAGGGGCAGTAAGCCCGTTGGGGATGGAGCCCAGCCCATGGCTGCCATGGGAGGCCTGAAGGTGCTTCTGCACTGGGCTGGTCCAGGCGGCGGGGAGCCCTGGGTCACTTTCAGTGAGTCATCGCTGACAGCTGAGGAAGTCTGCATCCACATTGCACATAAAGTTGGTATCACTCCTCCTTGCTTCAATCTCTTTGCCCTCTTCGATGCTCAGGCCCAAGTCTGGTTGCCCCCAAACCACATCCTAGAGATCCCCAGAGATGCAAGCCTGATGCTATATTTCCGCATAAGGTTTTATTTCCGGAACTGGCATGGCATGAATCCTCGGGAACCGGCTGTGTACCGTTGTGGGCCCCCAGGAACCGAGGCATCCTCAGATCAGACAGCACAGGGGATGCAACTCCTGGACCCAGCCTCATTTGAGTACCTCTTTGAGCAGGGCAAGCATGAGTTTGTGAATGACGTGGCATCACTGTGGGAGCTGTCGACCGAGGAGGAGATCCACCACTTTAAGAATGAGAGCCTGGGCATGGCCTTTCTGCACCTCTGTCACCTCGCTCTCCGCCATGGCATCCCCCTGGAGGAGGTGGCCAAGAAGACCAGCTTCAAGGACTGCATCCCGCGCTCCTTCCGCCGGCATATCCGGCAGCACAGCGCCCTGACCCGGCTGCGCCTTCGGAACGTCTTCCGCAGGTTCCTGCGGGACTTCCAGCCGGGCCGACTCTCCCAGCAGATGGTCATGGTCAAATACCTAGCCACACTCGAGCGGCTGGCACCCCGCTTCGGCACAGAGCGTGTGCCCGTGTGCCACCTGAGGCTGCTGGCCCAGGCCGAGGGGGAGCCCTGCTACATCCGGGACAGTGGGGTGGCCCCTACAGACCCTGGCCCTGAGTCTGCTGCTGGGCCCCCAACCCACGAGGTGCTGGTGACAGGCACTGGTGGCATCCAGTGGTGGCCAGTAGAGGAGGAGGTGAACAAGGAGGAGGGTTCTAGTGGCAGCAGTGGCAGGAACCCCCAAGCCAGCCTGTTTGGGAAGAAGGCCAAGGCTCACAAGGCAGTCGGCCAGCCGGCAGACAGGCCGCGGGAGCCACTGTGGGCCTACTTCTGTGACTTCCGGGACATCACCCACGTGGTGCTGAAAGAGCACTGTGTCAGCATCCACCGGCAGGACAACAAGTGCCTGGAGCTGAGCTTGCCTTCCCGGGCTGCGGCGCTGTCCTTCGTGTCGCTGGTGGACGGCTATTTCCGCCTGACGGCCGACTCCAGCCACTACCTGTGCCACGAGGTGGCTCCCCCACGGCTGGTGATGAGCATCCGGGATGGGATCCACGGACCCCTGCTGGAGCCATTTGTGCAGGCCAAGCTGCGGCCCGAGGACGGCCTGTACCTCATTCACTGGAGCACCAGCCACCCCTACCGCCTGATCCTCACAGTGGCCCAGCGTAGCCAGGCACCAGACGGCATGCAGAGCTTGCGGCTCCGAAAGTTCCCCATTGAGCAGCAGGACGGGGCCTTCGTGCTGGAGGGCTGGGGCCGGTCCTTCCCCAGCGTTCGGGAACTTGGGGCTGCCTTGCAGGGCTGCTTGCTGAGGGCCGGGGATGACTGCTTCTCTCTGCGTCGCTGTTGCCTGCCCCAACCAGGAGAAACCTCCAATCTCATCATCATGCGGGGGGCTCGGGCCAGCCCCAGGACACTCAACCTCAGCCAGCTCAGCTTCCACCGGGTTGACCAGAAGGAGATCACCCAGCTGTCCCACTTGGGCCAGGGCACAAGGACCAACGTGTATGAGGGCCGCCTGCGAGTGGAGGGCAGCGGGGACCCTGAGGAGGGCAAGATGGATGACGAGGACCCCCTCGTGCCTGGCAGGGACCGTGGGCAGGAGCTACGAGTGGTGCTCAAAGTGCTGGACCCTAGTCACCATGACATCGCCCTGGCCTTCTACGAGACAGCCAGCCTCATGAGCCAGGTCTCCCACACGCACCTGGCCTTCGTGCATGGCGTCTGTGTGCGCGGCCCTGAAAATATCATGGTGACAGAGTACGTGGAGCACGGACCCCTGGATGTGTGGCTGCGGAGGGAGCGGGGCCATGTGCCCATGGCTTGGAAGATGGTGGTGGCCCAGCAGCTGGCCAGCGCCCTCAGCTACCTGGAGAACAAGAACCTGGTTCATGGTAATGTGTGTGGCCGGAACATCCTGCTGGCCCGGCTGGGGTTGGCAGAGGGCACCAGCCCCTTCATCAAGCTGAGTGATCCTGGCGTGGGCCTGGGCGCCCTCTCCAGGGAGGAGCGGGTGGAGAGGATCCCCTGGCTGGCCCCCGAATGCCTACCAGGTGGGGCCAACAGCCTAAGCACCGCCATGGACAAGTGGGGGTTTGGCGCCACCCTCCTGGAGATCTGCTTTGACGGAGAGGCCCCTCTGCAGAGCCGCAGTCCCTCCGAGAAGGAGCATTTCTACCAGAGGCAGCACCGGCTGCCCGAGCCCTCCTGCCCACAGCTGGCCACACTCACCAGCCAGTGTCTGACCTATGAGCCAACCCAGAGGCCATCATTCCGCACCATCCTGCGTGACCTCACCCGGCTGCAGCCCCACAATCTTGCTGACGTCTTGACTGTGAACCCGGACTCACCGGCGTCGGACCCTACGGTTTTCCACAAGCGCTATTTGAAAAAGATCCGAGATCTGGGCGAGGGTCACTTCGGCAAGGTCAGCTTGTACTGCTACGATCCGACCAACGACGGCACTGGCGAGATGGTGGCGGTGAAAGCCCTCAAGGCAGACTGCGGCCCCCAGCACCGCTCGGGCTGGAAGCAGGAGATTGACATTCTGCGCACGCTCTACCACGAGCACATCATCAAGTACAAGGGCTGCTGCGAGGACCAAGGCGAGAAGTCGCTGCAGCTGGTCATGGAGTACGTGCCCCTGGGCAGCCTCCGAGACTACCTGCCCCGGCACAGCATCGGGCTGGCCCAGCTGCTGCTCTTCGCCCAGCAGATCTGCGAGGGCATGGCCTATCTGCACGCGCAGCACTACATCCACCGAGACCTAGCCGCGCGCAACGTGCTGCTGGACAACGACAGGCTGGTCAAGATCGGGGACTTTGGCCTAGCCAAGGCCGTGCCCGAAGGCCACGAGTACTACCGCGTGCGCGAGGATGGGGACAGCCCCGTGTTCTGGTATGCCCCAGAGTGCCTGAAGGAGTATAAGTTCTACTATGCGTCAGATGTCTGGTCCTTCGGGGTGACCCTGTATGAGCTGCTGACGCACTGTGACTCCAGCCAGAGCCCCCCCACGAAATTCCTTGAGCTCATAGGCATTGCTCAGGGTCAGATGACAGTTCTGAGACTCACTGAGTTGCTGGAACGAGGGGAGAGGCTGCCACGGCCCGACAAATGTCCCTGTGAGGTCTATCATCTCATGAAGAACTGCTGGGAGACAGAGGCGTCCTTTCGCCCAACCTTCGAGAACCTCATACCCATTCTGAAGACAGTCCATGAGAAGTACCAAGGCCAGGCCCCTTCAGTGTTCAGCGTGTGCTGA"
        ),
        "mutations": [
            ("Val15Ala", 43,45), ("Ala53Thr", 157,159), ("Arg197His", 589,591), ("Gly363Ser", 1087,1089),
            ("Arg744Trp", 2230,2232),
            ("Pro820His", 2458, 2460), ("His993Tyr", 2977, 2979), ("Glu1163Gly", 3487, 3489),
        ],
    },
}

# Fetch the standard codon table from Biopython
standard_table = CodonTable.unambiguous_dna_by_id[1]

# Define one-letter to three-letter conversion dictionary
aa_three_letter = {
    'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys',
    'E': 'Glu', 'Q': 'Gln', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
    'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro',
    'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val',
    '*': 'STOP'
}

# Convert Biopython's one-letter codon table to a three-letter format
codon_table = {codon: aa_three_letter[aa] for codon, aa in standard_table.forward_table.items()}

# Add stop codons
for stop_codon in standard_table.stop_codons:
    codon_table[stop_codon] = "STOP"


# Function to extract position from mutation string (e.g., "Val266Ile" → 266)
def extract_position(mutation_str):
    match = re.search(r'\d+', mutation_str)
    return int(match.group()) if match else None


# Function to translate DNA sequence into a protein sequence
def translate_dna(sequence):
    if len(sequence) < 3:
        return []

    translated_sequence = []  # List to store translated amino acids

    # Iterate through the sequence in steps of 3 (codons)
    for i in range(0, len(sequence), 3):
        codon = sequence[i:i + 3]  # Extract the triplet

        # Ensure it's a full triplet (skip last incomplete codons)
        if len(codon) == 3:
            amino_acid = codon_table.get(codon, 'Xxx')  # Translate to 3-letter code
            translated_sequence.append(amino_acid)

    return translated_sequence  # Return the translated amino acid sequence

# Extract surrounding nucleotide and translated amino acid sequences
flanking_sequences = []

with open(output_file, "w", newline="") as file:
    writer = csv.writer(file, delimiter="\t")
    writer.writerow(["Gene", "Mutation", "Upstream (Nuc)", "Mutated (Nuc)", "Downstream (Nuc)",
                     "Upstream (AA)", "Mutated (AA)", "Downstream (AA)", "Flanking Sequence (AA)"])

    for gene, data in gene_data.items():
        dna_sequence = data["sequence"]  # Full nucleotide sequence

        for mutation, start, end in data["mutations"]:
            # Ensure valid position within sequence bounds

            if start < 9 or end + 9 >= len(dna_sequence):
                continue  # Skip invalid cases

            if gene=="STAT1":
                # Extract 9 nucleotides upstream, 3 for mutation, 9 downstream
                mutated_nuc = dna_sequence[start: end+1]
                upstream_nuc = dna_sequence[start - 9: start]
                downstream_nuc = dna_sequence[end+1: end + 10]
            if gene!="STAT1":
                # Extract 9 nucleotides upstream, 3 for mutation, 9 downstream
                mutated_nuc = dna_sequence[start-1: end]
                upstream_nuc = dna_sequence[start - 10: start-1]
                downstream_nuc = dna_sequence[end: end + 9]

            # Translate nucleotide sequences to amino acids
            upstream_aa = translate_dna(upstream_nuc)
            mutated_aa = translate_dna(mutated_nuc)[0]  # Only 1 amino acid
            downstream_aa = translate_dna(downstream_nuc)

            # Format final flanking sequence
            flanking_seq = f"{'-'.join(upstream_aa)} [{mutated_aa}] {'-'.join(downstream_aa)}"

            # Store and write results
            flanking_sequences.append(flanking_seq)
            writer.writerow([gene, mutation, upstream_nuc, mutated_nuc, downstream_nuc,
                             '-'.join(upstream_aa), mutated_aa, '-'.join(downstream_aa), flanking_seq])

# Count the most common flanking sequences
flanking_counter = Counter(flanking_sequences)
most_common_flanking = flanking_counter.most_common(10)

# Print most common flanking sequences
print("Most common upstream and downstream sequences around mutations:")
for seq, count in most_common_flanking:
    print(f"{seq}: {count} occurrences")

print(f"Output saved to: {output_file}")
