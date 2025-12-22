"""
This file contains utilities for for phylogenetic tree analysis and modification

Author: Alex Poyer
"""

import subprocess
import re
from Bio import Phylo


def run_iqtree(input_alignment, output_prefix, iqtree_path="iqtree", extra_args=None):
    """
    Runs IQ-TREE on the given alignment file.

    Args:
        input_alignment (str): Path to the input alignment file (FASTA, PHYLIP, etc.).
        output_prefix (str): Prefix for output files.
        iqtree_path (str): Path to the IQ-TREE executable (default: "iqtree2").
        extra_args (list): Additional command-line arguments for IQ-TREE.

    Returns:
        int: The return code from the IQ-TREE process.
    """
    cmd = [iqtree_path, "-s", input_alignment, "-pre", output_prefix]
    if extra_args:
        cmd.extend(extra_args)
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
