import requests
import csv

# Input and output files
input_file = "cleaned_protein_ids.txt"
output_file = "taxonomy_info.csv"

# UniProt API base URL
uniprot_base_url = "https://rest.uniprot.org/uniprotkb/"

# Read protein IDs from the input file
with open(input_file, "r") as file:
    protein_ids = [line.strip() for line in file.readlines()]

taxonomy_data = []

# Fetch taxonomy information for each protein ID
for protein_id in protein_ids:
    try:
        # Send GET request to the UniProt API
        response = requests.get(uniprot_base_url + protein_id + ".json")
        response.raise_for_status()  # Raise an error for HTTP issues
        data = response.json()

        # Extract relevant taxonomy details
        taxonomy = data.get("organism", {})
        scientific_name = taxonomy.get("scientificName", "N/A")
        lineage = taxonomy.get("lineage", [])
        taxonomy_data.append([protein_id, scientific_name, " > ".join(lineage)])

        print(f"Processed: {protein_id}")
    except Exception as e:
        print(f"Error processing {protein_id}: {e}")
        taxonomy_data.append([protein_id, "Error", ""])

# Write taxonomy data to a CSV file
with open(output_file, "w", newline="") as file:
    writer = csv.writer(file)
    writer.writerow(["Protein ID", "Scientific Name", "Lineage"])
    writer.writerows(taxonomy_data)

print(f"Taxonomy information saved to {output_file}.")
