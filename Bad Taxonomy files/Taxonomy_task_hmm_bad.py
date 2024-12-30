import requests
import csv
import pandas as pd
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist
from collections import Counter
import numpy as np

# -----------------------------------------------------------------------------
# Step 1: Extract UniProt IDs from HMM Search Output
# -----------------------------------------------------------------------------

def extract_protein_ids_from_hmm(file_path):
    """
    Extract UniProt IDs from HMMER search output.
    Handles lines containing 'sp|' which indicate UniProt IDs.
    """
    protein_ids = []
    with open(file_path, "r") as file:
        for line in file:
            if "sp|" in line:  # UniProt IDs are marked with 'sp|'
                fields = line.strip().split()
                for field in fields:
                    if field.startswith("sp|"):  # Extract ID from 'sp|...' format
                        protein_id = field.split("|")[1]
                        protein_ids.append(protein_id)
    return list(set(protein_ids))  # Deduplicate IDs

# Specify HMM search output file and extract IDs
hmm_file = "hmmsearch_output.txt"
cleaned_protein_ids = extract_protein_ids_from_hmm(hmm_file)

# Save the cleaned IDs to a file
with open("hmm_cleaned_protein_ids.txt", "w") as output_file:
    output_file.write("\n".join(cleaned_protein_ids))

print(f"Extracted {len(cleaned_protein_ids)} protein IDs from {hmm_file}.")

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

output_file = "taxonomy_info_hmm.csv"
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
plt.title("Taxonomic Tree (HMM Search)", fontsize=16)
plt.xlabel("Distance", fontsize=14)
plt.ylabel("Lineage", fontsize=14)
plt.tight_layout()
plt.savefig("taxonomy_tree_hmm.png", dpi=300)
plt.show()

print("Taxonomic tree saved as taxonomy_tree_hmm.png.")

# -----------------------------------------------------------------------------
# Step 4: Taxonomic Tree with Node Sizes Proportional to Relative Abundance
# -----------------------------------------------------------------------------

# Calculate relative abundance
total_count = sum(lineage_counts.values())
relative_abundance = {k: v / total_count for k, v in lineage_counts.items()}

# Perform hierarchical clustering
condensed_distance_matrix = pdist(np.array(list(lineage_counts.values())).reshape(-1, 1), metric='euclidean')
linkage_matrix = linkage(condensed_distance_matrix, method='ward')

# Plot the dendrogram
plt.figure(figsize=(20, 40))
dendrogram_data = dendrogram(
    linkage_matrix,
    labels=lineages,
    orientation="right",
    leaf_font_size=10,
    color_threshold=0.5
)

# Add node sizes proportional to relative abundance
icoords = np.array(dendrogram_data['icoord'])
dcoords = np.array(dendrogram_data['dcoord'])
leaf_labels = np.array(dendrogram_data['ivl'])

# Map lineage to relative abundance for node sizes
lineage_to_abundance = {lineage: relative_abundance[lineage] * 1000 for lineage in leaf_labels}  # Scaled for visibility

# Plot scatter points for node sizes
for x, y, label in zip(icoords[:, 1], dcoords[:, 1], leaf_labels):
    if label in lineage_to_abundance:
        size = lineage_to_abundance[label]
        plt.scatter(x, y, s=size, color='blue', alpha=0.6, edgecolor='black', zorder=3)

plt.title("Taxonomic Tree with Node Sizes Proportional to Relative Abundance", fontsize=16)
plt.xlabel("Distance", fontsize=14)
plt.ylabel("Lineage", fontsize=14)
plt.tight_layout()

# Save and show plot
plt.savefig("taxonomy_tree_hmm_scaled.png", dpi=300)
plt.show()


# -----------------------------------------------------------------------------
# Step 5: Debugging Node Sizes
# -----------------------------------------------------------------------------

# Map lineage to relative abundance for node sizes
lineage_to_abundance = {lineage: relative_abundance[lineage] * 5000 for lineage in leaf_labels}  # Increased scaling factor

# Plot scatter points for node sizes
for x, y, label in zip(icoords[:, 1], dcoords[:, 1], leaf_labels):
    if label in lineage_to_abundance:
        size = lineage_to_abundance[label]
        plt.scatter(x, y, s=size, color='blue', alpha=0.7, edgecolor='black', zorder=3)

# Validate abundance values
print("Relative Abundance Values for Nodes:")
for lineage, abundance in relative_abundance.items():
    print(f"{lineage}: {abundance:.4f}")

plt.title("Taxonomic Tree with Node Sizes Proportional to Relative Abundance", fontsize=16)
plt.xlabel("Distance", fontsize=14)
plt.ylabel("Lineage", fontsize=14)
plt.tight_layout()

# Save and show plot
plt.savefig("taxonomy_tree_hmm_scaled_debug.png", dpi=300)
plt.show()


# Adjust scaling for node sizes
lineage_to_abundance = {lineage: relative_abundance[lineage] * 10000 for lineage in leaf_labels}  # Further scaled

# Plot scatter points for node sizes and validate
for x, y, label in zip(icoords[:, 1], dcoords[:, 1], leaf_labels):
    if label in lineage_to_abundance:
        size = lineage_to_abundance[label]
        plt.scatter(x, y, s=size, color='blue', alpha=0.6, edgecolor='black', zorder=3)
        # Add annotations to verify node size corresponds to abundance
        plt.text(x, y, f"{size:.2f}", fontsize=8, color='red', alpha=0.75)

# Validate abundance values
print("Relative Abundance Values for Nodes:")
for lineage, abundance in relative_abundance.items():
    print(f"{lineage}: {abundance:.4f}")
    
# Finalize and save the plot
plt.title("Taxonomic Tree with Node Sizes Proportional to Relative Abundance", fontsize=16)
plt.xlabel("Distance", fontsize=14)
plt.ylabel("Lineage", fontsize=14)
plt.tight_layout()
plt.savefig("taxonomy_tree_hmm_scaled_debug_2.png", dpi=300)
plt.show()