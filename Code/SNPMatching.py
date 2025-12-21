import sys
import pandas as pd
import re

def main(ped_file, recombinant_file, output_excel="matches_output.xlsx", ped_header="Position"):
    # Read ped.txt (space-separated)
    sample_data = {}
    positions = []
    with open(ped_file, 'r') as f:
        lines = f.readlines()
        if not lines:
            print("Error: ped.txt is empty")
            sys.exit(1)
        
        # Split first line by whitespace
        first_line = lines[0].strip()
        parts = re.split(r'\s+', first_line)
        parts = [p.strip() for p in parts]
        
        # Check header
        if parts[0].lower() == ped_header.lower():
            try:
                positions = [int(pos) for pos in parts[1:]]
                print(f"ped.txt positions (first 10): {positions[:10]}{'...' if len(positions) > 10 else ''}")
                print(f"Total positions in ped.txt: {len(positions)}")
            except ValueError:
                print(f"Invalid position values in ped.txt header: {parts[1:]}")
                sys.exit(1)
        else:
            print(f"Invalid ped.txt format: First column must be '{ped_header}', found '{parts[0]}'")
            sys.exit(1)
        
        # Subsequent lines: SampleX followed by 0/1 values
        for line in lines[1:]:
            parts = re.split(r'\s+', line.strip())
            parts = [p.strip() for p in parts]
            if len(parts) != len(positions) + 1:
                print(f"Skipping invalid line in ped.txt: {line.strip()}")
                continue
            sample_name = parts[0]
            try:
                values = [int(val) for val in parts[1:]]
                sample_data[sample_name] = dict(zip(positions, values))
                # Debug: Print SNPs for each sample
                snps = sorted([pos for pos, val in sample_data[sample_name].items() if val == 1])
                print(f"{sample_name} SNPs: {snps}")
                # Special check for Sample163
                if sample_name == 'Sample163':
                    expected_snps = [101, 109, 118, 139, 142, 205, 265, 277, 304, 334]
                    missing_snps = [snp for snp in expected_snps if snp not in snps]
                    extra_snps = [snp for snp in snps if snp not in expected_snps]
                    missing_from_ped = [snp for snp in expected_snps if snp not in positions]
                    print(f"Sample163: Missing expected SNPs: {missing_snps}")
                    print(f"Sample163: Extra unexpected SNPs: {extra_snps}")
                    print(f"Sample163: Expected SNPs not in ped.txt header: {missing_from_ped}")
                    sample163_positions_in_ped = [snp for snp in expected_snps if snp in positions]
                    print(f"Sample163: Expected SNPs in ped.txt header: {sample163_positions_in_ped}")
            except ValueError:
                print(f"Skipping line with non-integer values in ped.txt: {line.strip()}")
                continue

    # Read RecombinantSNPs.txt (space-separated)
    pairs = []
    all_l1_l2_positions = set()
    with open(recombinant_file, 'r') as f:
        lines = f.readlines()
        if not lines:
            print("Error: RecombinantSNPs.txt is empty")
            sys.exit(1)
        
        # Split header by whitespace
        first_line = lines[0].strip()
        header = re.split(r'\s+', first_line)
        header = [h.strip() for h in header]
        if header != ['L1', 'L2']:
            print(f"Invalid RecombinantSNPs.txt format: Header must be 'L1 L2', found {header}")
            sys.exit(1)
        
        # Process L1/L2 pairs
        for line in lines[1:]:
            parts = re.split(r'\s+', line.strip())
            parts = [p.strip() for p in parts]
            if len(parts) == 2:
                try:
                    l1 = int(float(parts[0]))
                    l2 = int(float(parts[1]))
                    pairs.append((l1, l2))
                    all_l1_l2_positions.add(l1)
                    all_l1_l2_positions.add(l2)
                except ValueError:
                    print(f"Skipping invalid line in RecombinantSNPs.txt: {line.strip()}")
                    continue
            else:
                print(f"Skipping line in RecombinantSNPs.txt with unexpected number of fields ({len(parts)}): {line.strip()}")
        # Debug: Print first few pairs and overlap with ped.txt
        print(f"First 5 L1/L2 pairs (up to {min(5, len(pairs))}): {pairs[:5]}")
        print(f"Total L1/L2 pairs: {len(pairs)}")
        # Check overlap between ped.txt positions and RecombinantSNPs.txt positions
        ped_positions = set(positions)
        overlap = ped_positions.intersection(all_l1_l2_positions)
        print(f"Positions in both ped.txt and RecombinantSNPs.txt: {sorted(overlap)}")
        if not overlap:
            print("Warning: No positions overlap between ped.txt and RecombinantSNPs.txt. No matches possible.")
        # Check Sample163 SNPs
        sample163_snps = [969, 991, 4633, 45396, 72601, 84847, 152294, 152391, 172834, 177693]
        pairs_with_163_snps = [(l1, l2) for l1, l2 in pairs if l1 in sample163_snps and l2 in sample163_snps]
        print(f"Pairs where both L1 and L2 are in Sample163 SNPs: {pairs_with_163_snps}")
        if not pairs_with_163_snps:
            print("Warning: No L1/L2 pairs in RecombinantSNPs.txt have both positions in Sample163's SNPs.")
        # Check potential matches for Sample163
        sample163_positions_in_ped = [p for p in sample163_snps if p in positions]
        potential_163_pairs = [(l1, l2) for l1, l2 in pairs if l1 in sample163_positions_in_ped and l2 in sample163_positions_in_ped]
        print(f"Potential Sample163 pairs (both L1 and L2 in ped.txt and Sample163 SNPs): {potential_163_pairs}")



    # Prepare data for output and comparison
    excel_data = []
    sample_to_pairs = {}
    for sample_name in sorted(sample_data.keys()):
        if not sample_name.startswith('Sample'):
            continue
        try:
            sample_num = int(sample_name[6:])
            if sample_num < 1 or sample_num > 213:
                continue
        except ValueError:
            continue  # Skip if not numeric after Sample
        
        sample_dict = sample_data[sample_name]
        matching_pairs = []
        # Check if both L1 and L2 have 1 in the sample
        for l1, l2 in pairs:
            if sample_dict.get(l1, 0) == 1 and sample_dict.get(l2, 0) == 1:
                matching_pairs.append((l1, l2))
        sample_to_pairs[sample_name] = tuple(sorted(matching_pairs))
        
    pairs_to_samples = {}
    for sample_name, matching_pairs in sample_to_pairs.items():
        pairs_to_samples.setdefault(matching_pairs, []).append(sample_name)
    
    for pair_set, samples in pairs_to_samples.items():
        num_pairs = len(pair_set)
        pair_identities = "; ".join([f"L1={l1}, L2={l2}" for l1, l2 in pair_set])
        unique = "Yes" if len(samples) == 1 else "No"
        shared_with = "" if len(samples) == 1 else ", ".join(samples)
        for sample in samples:
            excel_data.append([sample, num_pairs, unique, pair_identities, shared_with])
    # print(f"For {sample_name}:")
    # if matching_pairs:
        # for pair in matching_pairs:
            # print(f"  Match: L1={pair[0]}, L2={pair[1]}")
            # excel_data.append([sample_name, pair[0], pair[1]])
    # else:
        # print("  No matches found.")
        # excel_data.append([sample_name, "No matches found", ""])
    # print()

    # Save to Excel
    df = pd.DataFrame(excel_data, columns=["Sample", "Num Pairs", "UniqueRecombinant", "Pair Identities", "Shared With"])
    df.to_excel(output_excel, index=False)
    print(f"Results saved to {output_excel}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python SNPMatching.py ped.txt RecombinantSNPs.txt")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])