"""
This file contains the main workflow for the matching and phlyogenetic tree processes

Author: Alex Poyer
"""

import sys
import matching
import tree_utils
from pathlib import Path

ROOT_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT_DIR / "Data"
INPUT_DIR = DATA_DIR / "InputFiles"
OUTPUT_DIR = DATA_DIR / "OutputFiles"
IQTREE_OUT_DIR = DATA_DIR / "IQTree_out"
OUTPUT_PREFIX = "run_matching"
CONTREE = IQTREE_OUT_DIR / (OUTPUT_PREFIX + ".contree")

def main(
    ped_file,
    recombinant_file,
    alignment_file,
    output_excel=OUTPUT_DIR / "matches_output.xlsx",
):
    # tree_utils functions
    output_prefix = IQTREE_OUT_DIR / OUTPUT_PREFIX
    extra_args = ["-m", "GTR", "-bb", "1000", "-alrt", "1000", "-nt", "AUTO"]
    tree_utils.run_iqtree(alignment_file, output_prefix, extra_args=extra_args)
    time_pairs = tree_utils.parse_lengths(CONTREE)

    # matching functions
    sample_data = matching.read_ped_file(ped_file)
    pairs = matching.read_recombinant_file(recombinant_file)
    sample_to_pairs = matching.find_matching_pairs(sample_data, pairs)
    dataframe = matching.summarize_matches(sample_to_pairs, time_pairs)
    matching.save_to_excel(dataframe, output_excel)
    tree_utils.color_tree(dataframe, CONTREE, OUTPUT_DIR / (OUTPUT_PREFIX + ".nwk"))

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(
            "Usage: python run_matching.py <ped_file> <recombinant_file> <alignment_file>"
        )
        sys.exit(1)
    ped_file = sys.argv[1]
    recombinant_file = sys.argv[2]
    alignment_file = sys.argv[3]
    main(ped_file, recombinant_file, alignment_file)

# python run_matching.py ../Data/InputFiles/ped_12.22.2025.txt ../Data/InputFiles/HAPLOVIEWLDFILE ../Data/InputFiles/ALIGNEDSEQ213.fasta
