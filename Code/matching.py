"""
This file contains functions for reading PED and LD files and processing recombinant SNP pairs

Author: Alex Poyer
"""

import pandas as pd
import re
import tree_utils


# This function may further be optimized by utilizing pandas dataframes. Works for small datasets
def read_ped_file(ped_file, ped_header="Position"):
    """
    Reads a PED file and extracts genotype data for each sample, merging duplicate positions.

    Args:
        ped_file (str): Path to the PED text file.
        ped_header (str, optional): Expected header for the first column. Defaults to "Position".

    Raises:
        ValueError: If the PED file is empty or has an invalid format.

    Returns:
        dict: sample_data
            A dictionary mapping each sample name (str) to its genotype data.
            Duplicate positions are merged using logical OR.
            Format:
                {
                    "sample_name1": {position1: value1, position2: value2, ...},
                    "sample_name2": {position1: value1, position2: value2, ...},
                    ...
                }
            Where positions are SNP positions (int) and values are genotype calls (int, e.g., 0 or 1).
    """
    sample_data = {}
    with open(ped_file, "r") as f:
        lines = f.readlines()
        if not lines:
            raise ValueError("ped.txt is empty")
        first_line = lines[0].strip()
        parts = re.split(r"\s+", first_line)
        parts = [p.strip() for p in parts]

        if parts[0].lower() != ped_header.lower():
            errorString = f"Invalid ped.txt format: First column must be '{ped_header}', found '{parts[0]}'"
            raise ValueError(errorString)

        # Detect duplicate positions
        positions = [int(pos) for pos in parts[1:]]
        pos_indices = {}
        for idx, pos in enumerate(positions):
            pos_indices.setdefault(pos, []).append(idx)

        for line in lines[1:]:
            parts = re.split(r"\s+", line.strip())
            parts = [p.strip() for p in parts]
            if len(parts) != len(positions) + 1:
                print(f"Skipping invalid line in ped.txt: {line.strip()}")
                continue
            sample_name = parts[0]
            try:
                values = [int(val) for val in parts[1:]]
                merged = {}
                for pos, idxs in pos_indices.items():
                    # Merge duplicate columns using logical OR
                    merged[pos] = int(any(values[i] for i in idxs))
                sample_data[sample_name] = merged
            except ValueError:
                continue
    return sample_data


# primative approach for parsing and storing data -> efficient for small/medium datasets
def read_recombinant_file(recombinant_file, CIHi_filter=None):
    """
    Parses a Haploview LD export file and returns all recombinant SNP position pairs.

    Args:
        recombinant_file (str): Path to the LD text file.
        CIHi_filter (float): minimum CIHi for a pair to be considered

    Raises:
        ValueError: If the file does not contain columns L1, L2, and CIhi.

    Returns:
        List: pairs
            A list of recombinant SNP position pairs to check for in each sample.
            Format:
                [
                    (l1_1, l2_1),
                    (l1_2, l2_2),
                    ...
                ]
            Where each tuple contains two SNP positions (int).
    """
    df = pd.read_csv(recombinant_file, sep=r"\t", engine="python", usecols=range(9))
    # print(df.columns) # debug
    if not {"L1", "L2", "CIhi"}.issubset(df.columns):
        raise ValueError("recombinant_file must contain columns: L1, L2, CIhi")
    pairs = []
    for _, row in df.iterrows():
        # print(row) # debug
        try:
            l1 = int(row["L1"])
            l2 = int(row["L2"])
            cihi = float(row["CIhi"])
            if CIHi_filter and cihi < CIHi_filter:
                continue
            pairs.append((l1, l2))
        except (ValueError, KeyError):
            continue
    return pairs


def find_matching_pairs(sample_data, pairs):
    """Finds which recombinant SNP pairs are present in each sample

    Args:
        sample_data (dict): A dictionary mapping each sample name (str) to its genotype data.
            See read_ped_file for format.
        pairs (list): A list of recombinant SNP position pairs to check for in each sample.
            See read_recombinant_file for format.

    Returns:
        dict: sample_to_pairs
            A dictionary mapping each sample name (str) to the set of recombinant pairs present in that sample.
            Format:
                {
                    "sample_name1": ((l1_1, l2_1), (l1_2, l2_2), ...),
                    "sample_name2": ((l1_3, l2_3), ...),
                    ...
                }
            Where each value is a tuple of recombinant pairs (each pair is a tuple of two SNP positions, int).
    """
    sample_to_pairs = {}
    for sample_name in sorted(sample_data.keys()):
        # print(f"matching {sample_name}") # Debug
        sample_dict = sample_data[sample_name]
        matching_pairs = []
        for l1, l2 in pairs:
            if sample_dict.get(l1, 0) == 1 and sample_dict.get(l2, 0) == 1:
                matching_pairs.append((l1, l2))
        sample_to_pairs[sample_name] = tuple(sorted(matching_pairs))
    # print("pairs:", sample_to_pairs) # Debug
    return sample_to_pairs


def summarize_matches(sample_to_pairs, phylo_times=None):
    """
    Summarizes which samples share the same set of recombinant pairs.

    Args:
        sample_to_pairs (dict): A dictionary mapping each sample name (str) to the set of recombinant pairs present in that sample.
            See find_matching_pairs for format.
        phylo_times (dict, optional): A dictionary mapping each sample name (str) to branch length / phylogenetic time (int)

    Returns:
        DataFrame: df
            A dataframe containing rows :
                [epi_isl identifier, sample name, number of pairs, uniqueness, pair identities, samples it shares identities with]
    """
    pairs_to_samples = {}
    for sample_name, matching_pairs in sample_to_pairs.items():
        pairs_to_samples.setdefault(matching_pairs, []).append(sample_name)
    rows = []
    for pair_set, samples in pairs_to_samples.items():
        num_pairs = len(pair_set)
        pair_identities = "; ".join([f"L1={l1}, L2={l2}" for l1, l2 in pair_set])
        unique = "Yes" if len(samples) == 1 else "No"
        shared_with = "" if len(samples) == 1 else ", ".join(samples)
        for sample in samples:
            epi_isl = tree_utils.extract_epi_isl(sample)
            phylo_time = phylo_times.get(epi_isl, "") if phylo_times else ""
            rows.append(
                [
                    epi_isl,
                    sample,
                    num_pairs,
                    unique,
                    pair_identities,
                    shared_with,
                    phylo_time,
                ]
            )
    df = pd.DataFrame(
        rows,
        columns=[
            "epi_isl",
            "Sample",
            "Num Pairs",
            "UniqueRecombinant",
            "Pair Identities",
            "Shared With",
            "Phylogenetic Time",
        ],
    )
    return df


def save_to_excel(df, output_excel):
    """
    Saves the DataFrame to an Excel file


    Args:
        excel_data (list): list with summarized sample data
        output_excel (str): output excel file name
    """
    export_cols = [
        "Sample",
        "Num Pairs",
        "UniqueRecombinant",
        "Pair Identities",
        "Shared With",
        "Phylogenetic Time",
    ]
    df.to_excel(output_excel, index=False, columns=export_cols)
    print(f"Results saved to {output_excel}")
