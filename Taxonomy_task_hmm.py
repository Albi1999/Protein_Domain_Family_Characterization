import requests
import csv
import pandas as pd
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist
from collections import Counter
import numpy as np
from typing import List, Dict, Tuple
import time
from tqdm import tqdm
import logging
import os
from datetime import datetime

class TaxonomyAnalyzer:
    def __init__(self, hmm_file: str, max_retries: int = 3, retry_delay: int = 1):
        self.hmm_file = hmm_file
        self.max_retries = max_retries
        self.retry_delay = retry_delay
        self.uniprot_base_url = "https://rest.uniprot.org/uniprotkb/"
        
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            filename="Final_taxonomy_analysis.log"
        )

    def extract_protein_ids(self) -> List[str]:
        """
        Extract UniProt IDs from HMMER search output with improved error handling.
        """
        try:
            protein_ids = []
            with open(self.hmm_file, "r") as file:
                for line in file:
                    if "sp|" in line:
                        fields = line.strip().split()
                        for field in fields:
                            if field.startswith("sp|"):
                                protein_id = field.split("|")[1]
                                protein_ids.append(protein_id)

            unique_ids = list(set(protein_ids))
            logging.info(f"Extracted {len(unique_ids)} unique protein IDs")
            return unique_ids
        
        except FileNotFoundError:
            logging.error(f"HMM file not found: {self.hmm_file}")
            raise
        except Exception as e:
            logging.error(f"Error extracting protein IDs: {str(e)}")
            raise

    def fetch_taxonomy_info(self, protein_ids: List[str], output_file: str) -> str:
        """
        Fetch taxonomy information with improved error handling and progress tracking.
        """
        taxonomy_data = []
        error_counts = {"success": 0, "failed": 0}

        for protein_id in tqdm(protein_ids, desc="Fetching taxonomy info"):
            for attempt in range(self.max_retries):
                try:
                    response = requests.get(f"{self.uniprot_base_url}{protein_id}.json")
                    response.raise_for_status()
                    data = response.json()

                    taxonomy = data.get("organism", {})
                    scientific_name = taxonomy.get("scientificName", "N/A")
                    lineage = taxonomy.get("lineage", [])
                    taxonomy_data.append([protein_id, scientific_name, " > ".join(lineage)])

                    error_counts["success"] += 1
                    break
                
                except requests.exceptions.RequestException as e:
                    if attempt == self.max_retries - 1:
                        logging.error(f"Failed to fetch {protein_id} after {self.max_retries} attempts: {str(e)}")
                        taxonomy_data.append([protein_id, "Error", ""])
                        error_counts["failed"] += 1
                    else:
                        time.sleep(self.retry_delay)

        with open(output_file, "w", newline="") as file:
            writer = csv.writer(file)
            writer.writerow(["Protein ID", "Scientific Name", "Lineage"])
            writer.writerows(taxonomy_data)
        logging.info(f"Taxonomy data saved to {output_file}")

        return output_file

    def create_taxonomic_tree(self, taxonomy_file: str, output_file: str, 
                              lineage_depth: int = 3, figsize: Tuple[int, int] = (15, 20),
                              min_node_size: int = 50, max_node_size: int = 500) -> str:
        """
        Create an improved taxonomic tree visualization.
        """
        # Load and process taxonomy data
        taxonomy_info = pd.read_csv(taxonomy_file)
        filtered_taxonomy = taxonomy_info[taxonomy_info['Scientific Name'] != "Error"]

        filtered_taxonomy['Simplified_Lineage'] = filtered_taxonomy['Lineage'].apply(
            lambda x: " > ".join(x.split(" > ")[:lineage_depth])
        )

        # Calculate abundance with improved scaling
        lineage_counts = Counter(filtered_taxonomy['Simplified_Lineage'])
        lineage_df = pd.DataFrame(list(lineage_counts.items()), columns=['Lineage', 'Count'])

        max_count = lineage_df['Count'].max()
        lineage_df['Scaled_Size'] = (lineage_df['Count'] / max_count) * (max_node_size - min_node_size) + min_node_size

        # Perform hierarchical clustering
        lineages = lineage_df['Lineage'].tolist()
        counts = lineage_df['Count'].values.reshape(-1, 1)

        distance_matrix = pdist(counts, metric='euclidean')
        linkage_matrix = linkage(distance_matrix, method='ward')

        # Create plot
        fig, ax = plt.subplots(figsize=figsize)

        dendrogram_data = dendrogram(
            linkage_matrix,
            labels=lineages,
            orientation="right",
            leaf_font_size=10,
            leaf_rotation=0,
            color_threshold=0.7 * max(linkage_matrix[:, 2]),
            distance_sort='descending',
            ax=ax
        )

        norm = plt.Normalize(vmin=min(lineage_df['Count']), vmax=max(lineage_df['Count']))

        # Add nodes
        for x, y, label in zip(dendrogram_data['icoord'], dendrogram_data['dcoord'], dendrogram_data['ivl']):
            if label in lineage_df['Lineage'].values:
                count = lineage_df.loc[lineage_df['Lineage'] == label, 'Count'].iloc[0]
                size = lineage_df.loc[lineage_df['Lineage'] == label, 'Scaled_Size'].iloc[0]
                color = plt.cm.viridis(norm(count))
                ax.scatter(x[1], y[1], s=size, color=color, edgecolor='black', zorder=3)

        ax.set_title("Taxonomic Tree with Abundance Visualization", fontsize=16)
        ax.set_xlabel("Hierarchical Distance", fontsize=12)
        ax.set_ylabel("Taxonomic Lineage", fontsize=12)

        # Add colorbar
        sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, pad=0.02)
        cbar.set_label('Species Count', fontsize=12)

        plt.tight_layout()
        plt.savefig(output_file, dpi=300)
        logging.info(f"Taxonomic tree saved to {output_file}")

        return output_file

def main():
    """Main execution function."""
    analyzer = TaxonomyAnalyzer("hmmsearch_output.txt")

    protein_ids = analyzer.extract_protein_ids()
    taxonomy_file = analyzer.fetch_taxonomy_info(protein_ids, "Final_taxonomy_info.csv")
    tree_file = analyzer.create_taxonomic_tree(
        taxonomy_file, "Final_taxonomy_tree.png", lineage_depth=3, figsize=(15, 8)
    )

    print("\nFiles created:")
    print(f"Taxonomy data: {taxonomy_file}")
    print(f"Taxonomic tree: {tree_file}")

if __name__ == "__main__":
    main()
