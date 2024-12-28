import pandas as pd
from collections import Counter
from scipy.cluster.hierarchy import dendrogram, linkage
import matplotlib.pyplot as plt

# Load taxonomy data
taxonomy_file = "taxonomy_info.csv"  # Replace with the correct file path
taxonomy_data = pd.read_csv(taxonomy_file)

# Ensure the "Lineage" column is processed correctly
taxonomy_data['Lineage'] = taxonomy_data['Lineage'].fillna("Unknown")

# Count occurrences of each lineage
lineage_counts = Counter(taxonomy_data['Lineage'])

# Create a matrix for hierarchical clustering
lineage_labels = list(lineage_counts.keys())
lineage_sizes = [[count] for count in lineage_counts.values()]

# Perform hierarchical clustering
Z = linkage(lineage_sizes, method='ward')

# Adjust figure size (in inches) and resolution (DPI)
fig_width = 10  # Adjust width
fig_height = 6  # Adjust height
dpi = 300       # Set high DPI for better quality

# Plot the dendrogram
plt.figure(figsize=(fig_width, fig_height), dpi=dpi)
dendrogram(Z, labels=lineage_labels, orientation='right', leaf_font_size=8)
plt.title("Taxonomy Tree", fontsize=14)
plt.xlabel("Counts", fontsize=12)
plt.ylabel("Lineage", fontsize=12)
plt.tight_layout()

# Save the plot
output_file = "tree_output.png"  # Set your desired output file name
plt.savefig(output_file, dpi=dpi)
plt.show()
