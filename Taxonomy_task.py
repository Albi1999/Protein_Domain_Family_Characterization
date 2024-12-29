# Taxonomy Task: Extract UniProt IDs from a FASTA file, 
# fetch taxonomy information for each ID using the UniProt API, 
# and visualize the taxonomic tree based on the lineage information.

import requests
import csv
import pandas as pd
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist
from collections import Counter
import numpy as np

# -----------------------------------------------------------------------------
# Step 1: Extract UniProt IDs from a FASTA file
# -----------------------------------------------------------------------------

def extract_uniprot_ids(fasta_file):
    """Extract UniProt IDs from a FASTA file."""
    cleaned_ids = []
    with open(fasta_file, "r") as file:
        for line in file:
            if line.startswith(">"):
                protein_id = line.split()[0][1:]  # Remove '>' and take first token
                if "|" in protein_id:
                    protein_id = protein_id.split("|")[1]  # Extract ID between '|'
                cleaned_ids.append(protein_id)
    return cleaned_ids

# Specify FASTA file and extract IDs
fasta_file = "trimmed_alignment.fasta"
cleaned_protein_ids = extract_uniprot_ids(fasta_file)

# Save the cleaned IDs to a file
with open("cleaned_protein_ids.txt", "w") as output_file:
    output_file.write("\n".join(cleaned_protein_ids))

print(f"Extracted {len(cleaned_protein_ids)} UniProt IDs.")

# -----------------------------------------------------------------------------
# Step 2: Fetch Taxonomy Information for Protein IDs
# -----------------------------------------------------------------------------

def fetch_taxonomy_info(protein_ids, output_file):
    """Fetch taxonomy information for a list of UniProt IDs and save to CSV."""
    uniprot_base_url = "https://rest.uniprot.org/uniprotkb/"
    taxonomy_data = []

    for protein_id in protein_ids:
        try:
            response = requests.get(uniprot_base_url + protein_id + ".json")
            response.raise_for_status()
            data = response.json()

            taxonomy = data.get("organism", {})
            scientific_name = taxonomy.get("scientificName", "N/A")
            lineage = taxonomy.get("lineage", [])
            taxonomy_data.append([protein_id, scientific_name, " > ".join(lineage)])

            print(f"Processed: {protein_id}")
        except requests.exceptions.RequestException as e:
            print(f"Error processing {protein_id}: {e}")
            taxonomy_data.append([protein_id, "Error", ""])

    # Write taxonomy data to CSV
    with open(output_file, "w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["Protein ID", "Scientific Name", "Lineage"])
        writer.writerows(taxonomy_data)

output_file = "taxonomy_info.csv"
fetch_taxonomy_info(cleaned_protein_ids, output_file)
print(f"Taxonomy information saved to {output_file}.")

# -----------------------------------------------------------------------------
# Step 3: Taxonomic Tree Visualization
# -----------------------------------------------------------------------------

def simplify_lineage(lineage):
    """Simplify lineage to keep only genus and species."""
    levels = lineage.split(" > ")
    if len(levels) >= 2:
        return " > ".join(levels[-2:])
    return lineage

# Load taxonomy information
taxonomy_info = pd.read_csv(output_file)
filtered_taxonomy = taxonomy_info[taxonomy_info['Scientific Name'] != "Error"]
filtered_taxonomy['Simplified_Lineage'] = filtered_taxonomy['Lineage'].apply(simplify_lineage)

# Count the abundance of each taxonomic group
lineage_counts = Counter(filtered_taxonomy['Simplified_Lineage'])
lineage_df = pd.DataFrame(list(lineage_counts.items()), columns=['Lineage', 'Count'])

# Perform hierarchical clustering
lineages = lineage_df['Lineage'].tolist()
counts = lineage_df['Count'].tolist()
condensed_distance_matrix = pdist(np.array(counts).reshape(-1, 1), metric='euclidean')
linkage_matrix = linkage(condensed_distance_matrix, method='ward')

# Plot the dendrogram
plt.figure(figsize=(20, 40))
dendrogram(
    linkage_matrix,
    labels=lineages,
    orientation="right",
    leaf_font_size=10,
    color_threshold=0.5
)
plt.title("Taxonomic Tree", fontsize=16)
plt.xlabel("Distance", fontsize=14)
plt.ylabel("Lineage", fontsize=14)
plt.tight_layout()
plt.savefig("taxonomy_tree.png", dpi=300)
plt.show()

print("Taxonomic tree saved as taxonomy_tree.png.")
