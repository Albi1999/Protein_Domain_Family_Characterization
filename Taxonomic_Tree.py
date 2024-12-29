import pandas as pd
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist, squareform
from collections import Counter
import numpy as np

# Load cleaned protein IDs
with open("cleaned_protein_ids.txt", "r") as f:
    protein_ids = f.read().splitlines()

# Load taxonomy info
taxonomy_info = pd.read_csv("taxonomy_info.csv")

# Filter taxonomy data to include only relevant protein IDs
filtered_taxonomy = taxonomy_info[taxonomy_info['Protein ID'].isin(protein_ids)]

# Simplify taxonomic lineage (e.g., keeping genus and species)
def simplify_lineage(lineage):
    levels = lineage.split(" > ")
    if len(levels) >= 2:
        return " > ".join(levels[-2:])  # Keep genus and species
    return lineage

filtered_taxonomy['Simplified_Lineage'] = filtered_taxonomy['Lineage'].apply(simplify_lineage)

# Count the abundance of each taxonomic group
lineage_counts = Counter(filtered_taxonomy['Simplified_Lineage'])

# Convert counts to a DataFrame
lineage_df = pd.DataFrame(list(lineage_counts.items()), columns=['Lineage', 'Count'])

# Perform hierarchical clustering
lineages = lineage_df['Lineage'].tolist()
counts = lineage_df['Count'].tolist()

# Create a condensed distance matrix using pdist
condensed_distance_matrix = pdist(np.array(counts).reshape(-1, 1), metric='euclidean')
linkage_matrix = linkage(condensed_distance_matrix, method='ward')

# Plot the dendrogram
plt.figure(figsize=(15, 30))  # Adjust figure size for better readability

# Plot dendrogram
dendrogram(
    linkage_matrix,
    labels=lineages,
    orientation="right",
    leaf_font_size=8,  # Smaller font size for better spacing
    color_threshold=0.5
)

# Adjust plot aesthetics
plt.title("Taxonomy Tree with Scaled Node Sizes", fontsize=14)
plt.xlabel("Distance", fontsize=12)
plt.ylabel("Lineage", fontsize=12)
plt.tight_layout()

# Save and display the plot
plt.savefig("taxonomy_tree.png", dpi=300)
plt.show()
