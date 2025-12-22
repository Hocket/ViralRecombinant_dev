import subprocess


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
