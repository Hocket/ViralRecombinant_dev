"""
This file contains the main workflow for the matching and phlyogenetic tree processes

Author: Alex Poyer
"""

import sys
import matching
import tree_utils


def main(
    ped_file,
    recombinant_file,
    alignment_file,
    output_excel="../Data/OutputFiles/matches_output.xlsx",
):
    # tree_utils functions
    output_prefix = "../Data/IQTree_out/run_matching"
    extra_args = ["-m", "GTR", "-bb", "1000", "-alrt", "1000", "-nt", "AUTO"]
    tree_utils.run_iqtree(alignment_file, output_prefix, extra_args=extra_args)
    time_pairs = tree_utils.parse_lengths("../Data/IQTree_out/run_matching.treefile")

    # matching functions
    sample_data = matching.read_ped_file(ped_file)
    pairs = matching.read_recombinant_file(recombinant_file)
    sample_to_pairs = matching.find_matching_pairs(sample_data, pairs)
    excel_data = matching.summarize_matches(sample_to_pairs, time_pairs)
    matching.save_to_excel(excel_data, output_excel)


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
