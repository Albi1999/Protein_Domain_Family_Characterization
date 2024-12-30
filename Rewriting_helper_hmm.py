import re
import csv

# File paths
input_file_path = "hmmsearch_output.txt"
output_file_path = "hmmsearch_output.csv"

# Initialize storage for parsed data
parsed_data = []

# Regular expressions to capture key information
header_regex = r">> ([^\s]+)"
domain_regex = r"\s+(\d+) [!?]\s+([\d\.]+)\s+[\d\.]+\s+[\de\.\+\-]+\s+[\de\.\+\-]+\s+\d+\s+\d+\s+(?:\[\.|\.\.)+\s+(\d+)\s+(\d+)"

with open(input_file_path, "r") as infile:
    current_protein = None

    for line in infile:
        # Match protein header line
        header_match = re.match(header_regex, line)
        if header_match:
            # If we already captured a protein, save its data
            if current_protein:
                parsed_data.append(current_protein)

            # Start a new protein record
            protein_id = header_match.groups()[0]
            current_protein = {
                "protein_name": protein_id.split("|")[2],
                "uniprot_id": protein_id.split("|")[1],
                "domains": []
            }

        # Match domain annotation (including both `!` and `?` lines)
        domain_match = re.match(domain_regex, line)
        if domain_match and current_protein:
            _, score, start, end = domain_match.groups()
            start, end, score = int(start), int(end), float(score)
            length = end - start + 1
            current_protein["domains"].append((score, start, end, length))

    # Handle the last protein record
    if current_protein:
        parsed_data.append(current_protein)

# Prepare fieldnames dynamically
fieldnames = ["protein_name", "uniprot_id"]
max_domains = max(len(protein["domains"]) for protein in parsed_data)
for i in range(1, max_domains + 1):
    fieldnames.extend([
        f"domain_{i}_score", f"domain_{i}_start", f"domain_{i}_end", f"domain_{i}_length"
    ])

# Write to CSV
with open(output_file_path, "w", newline="") as outfile:
    writer = csv.DictWriter(outfile, fieldnames=fieldnames)
    writer.writeheader()
    for protein in parsed_data:
        row = {
            "protein_name": protein["protein_name"],
            "uniprot_id": protein["uniprot_id"]
        }
        for i, domain in enumerate(protein["domains"], start=1):
            row[f"domain_{i}_score"] = domain[0]
            row[f"domain_{i}_start"] = domain[1]
            row[f"domain_{i}_end"] = domain[2]
            row[f"domain_{i}_length"] = domain[3]
        writer.writerow(row)

print(f"CSV file generated: {output_file_path}")
