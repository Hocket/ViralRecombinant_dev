import os

# Input and output files
mapping_file = "character_mapping.txt"
query_file = "position.txt"
output_file = "position_values.txt"

# Check if input files exist
if not os.path.isfile(mapping_file):
    print(f"Error: {mapping_file} not found")
    exit(1)
if not os.path.isfile(query_file):
    print(f"Error: {query_file} not found")
    exit(1)

# Load character mapping
mapping = {}
with open(mapping_file, 'r') as f:
    for line in f:
        pos, val = line.strip().split(':')
        mapping[int(pos)] = val

# Process query positions and write to output file
with open(output_file, 'w') as out:
    with open(query_file, 'r') as f:
        for line in f:
            query_pos = line.strip()
            if not query_pos.isdigit():
                out.write("Invalid\n")
                continue
            query_pos = int(query_pos)

            if query_pos not in mapping:
                out.write("NotFound\n")
                continue

            value = mapping[query_pos]
            if value == '-':
                # Find flanking alphabetical characters
                prev_value = "None"
                next_value = "None"
                
                # Search backward
                for i in range(query_pos - 1, 0, -1):
                    if i in mapping and mapping[i] != '-':
                        prev_value = mapping[i]
                        break
                
                # Search forward
                for i in range(query_pos + 1, max(mapping.keys()) + 1):
                    if i in mapping and mapping[i] != '-':
                        next_value = mapping[i]
                        break
                
                # Write flanking values in number-number format, or NotFound if both are None
                if prev_value == "None" and next_value == "None":
                    out.write("NotFound\n")
                else:
                    out.write(f"{prev_value}-{next_value}\n")
            else:
                out.write(f"{value}\n")
