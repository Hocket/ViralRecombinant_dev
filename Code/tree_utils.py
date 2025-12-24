"""
This file contains utilities for for phylogenetic tree analysis and modification

Author: Alex Poyer
"""

import subprocess
import re
from Bio import Phylo


def run_iqtree(input_alignment, output_prefix, iqtree_path="iqtree3", extra_args=None):
    """
    Runs IQ-TREE on the given alignment file.

    Args:
        input_alignment (str): Path to the input alignment file (FASTA, PHYLIP, etc.).
        output_prefix (str): Prefix for output files.
        iqtree_path (str): Path to the IQ-TREE executable (default: "iqtree3").
        extra_args (list): Additional command-line arguments for IQ-TREE.

    Returns:
        int: The return code from the IQ-TREE process.
    """
    cmd = [str(iqtree_path), "-s", str(input_alignment), "-pre", str(output_prefix)]
    if extra_args:
        cmd.extend(map(str,extra_args))
    print("Running:", " ".join(cmd))
    result = subprocess.run(cmd)
    return result.returncode


def parse_lengths(tree_file):
    tree = Phylo.read(tree_file, "newick")
    time_pairs = {}
    for clade in tree.get_terminals():
        epi_isl = extract_epi_isl(clade.name)
        if epi_isl:
            time_pairs[epi_isl] = clade.branch_length
    return time_pairs


def extract_epi_isl(name):
    """Extracts the EPI_ISL identifier from a sample name string."""
    match = re.search(r"EPI_ISL_\d+", name)
    return match.group(0) if match else None


def color_tree(data, tree_file, output_file):
    """ writes comments on branches to state color

    Args:
        data (dataFrame): pandas dataframe representing a tree
        tree_file (string): A newick style tree file to be colored
        output_file (string): The output file for the colored newick tree
    """

    # pass through and add color to tree
    tree = Phylo.read(tree_file, "newick")
    
    df = data.set_index("epi_isl")
    for clade in tree.get_terminals():
        sample_name = clade.name
        epi_isl = extract_epi_isl(sample_name)
        if epi_isl is None:
            continue
        unique_recombinant = (
            df.at[epi_isl, "UniqueRecombinant"]
            if epi_isl in df.index
            else None
        )
        # COLORS USE HEX CODE
        if unique_recombinant.lower() == "yes":
            clade.comment = "&color=#FF0000"
            continue
        # TO IMPLEMENT WHEN COLORS ARE KNOWN
        # num_pairs = (
        #     df.at[epi_isl, "Num Pairs"]
        #     if epi_isl in df.index
        #     else None
        # )
    Phylo.write(tree, output_file, "newick")
    print("Color tree saved to", output_file)