import requests
import pandas as pd
from ete3 import Tree, TreeStyle, NodeStyle, TextFace
from tqdm import tqdm
import time

# TaxonomyAnalyzer Class for fetching taxonomy information
class TaxonomyAnalyzer:
    def __init__(self, max_retries: int = 3, retry_delay: int = 1):
        self.max_retries = max_retries
        self.retry_delay = retry_delay
        self.uniprot_base_url = "https://rest.uniprot.org/uniprotkb/"

    def fetch_taxonomy_info(self, protein_ids: list, output_file: str):
        taxonomy_data = []

        pbar = tqdm(protein_ids, desc="Fetching taxonomy data")

        for protein_id in pbar:
            pbar.set_description(f"Processing {protein_id}")

            for attempt in range(self.max_retries):
                try:
                    response = requests.get(f"{self.uniprot_base_url}{protein_id}.json")
                    response.raise_for_status()
                    data = response.json()

                    taxonomy = data.get("organism", {})
                    scientific_name = taxonomy.get("scientificName", "N/A")
                    lineage = taxonomy.get("lineage", [])

                    taxonomy_data.append([protein_id, scientific_name, " > ".join(lineage)])
                    break

                except requests.exceptions.RequestException as e:
                    print(f"Error fetching data for {protein_id}: {e}")
                    if attempt == self.max_retries - 1:
                        taxonomy_data.append([protein_id, "Error", ""])
                    else:
                        time.sleep(self.retry_delay)

        taxonomy_df = pd.DataFrame(taxonomy_data, columns=["Protein ID", "Scientific Name", "Lineage"])
        taxonomy_df.to_csv(output_file, index=False)
        return taxonomy_df

# Load protein IDs from files
def load_protein_ids(psiblast_file, hmm_file, e_threshold=0.001):
    psiblast_df = pd.read_csv(psiblast_file)
    hmm_df = pd.read_csv(hmm_file)

    filtered_hmm_proteins = hmm_df[hmm_df['E-value'] <= e_threshold]['uniprot_id']
    psiblast_proteins = set(psiblast_df['uniprot_id'])
    hmm_proteins = set(filtered_hmm_proteins)

    return list(psiblast_proteins.union(hmm_proteins))

# Fetch taxonomy data
def main():
    psiblast_file = "psiblast_parsed.csv"
    hmm_file = "hmmsearch_output.csv"
    protein_ids = load_protein_ids(psiblast_file, hmm_file)

    analyzer = TaxonomyAnalyzer()
    taxonomy_data = analyzer.fetch_taxonomy_info(protein_ids, "taxonomy_info.csv")

    # Check and print column names
    print("Columns in the file:", taxonomy_data.columns)

    # Correct column name
    correct_column_name = "Lineage"  # Use the correct column name

    # Process taxonomy data
    def process_taxonomy(data):
        taxonomy_dict = {}
        for _, row in data.iterrows():
            lineage = row[correct_column_name].split(" > ")  # Ensure semicolon-separated lineage
            current = taxonomy_dict
            for level in lineage:
                if level not in current:
                    current[level] = {}
                current = current[level]
        return taxonomy_dict

    # Create a nested dictionary of taxonomy
    taxonomy_dict = process_taxonomy(taxonomy_data)

    # Count relative abundance
    abundance_counts = taxonomy_data[correct_column_name].value_counts().to_dict()

    # Create a Newick string for the taxonomy tree
    def dict_to_newick(d, parent_abundance=None):
        newick = ""
        for key, sub_dict in d.items():
            size = parent_abundance.get(key, 1) if parent_abundance else 1
            sub_tree = dict_to_newick(sub_dict, parent_abundance)
            newick += f"({sub_tree}){key}:{size}," if sub_tree else f"{key}:{size},"
        return newick.rstrip(",")

    newick_tree = f"({dict_to_newick(taxonomy_dict, abundance_counts)});"

    # Plot using ETE Toolkit
    phylo_tree = Tree(newick_tree, format=1)
    tree_style = TreeStyle()
    tree_style.show_leaf_name = False

    # Adjust node labels to show lineage at all nodes
    for node in phylo_tree.traverse():
        node.add_face(TextFace(node.name, fsize=10, fgcolor="black"), column=0)

    # Adjust node sizes (normalize and refine scaling)
    max_size = 50  # Increase max size for better differentiation
    scaling_factor = 2  # Further refine scaling for visual contrast
    for node in phylo_tree.traverse():
        nstyle = NodeStyle()
        size = abundance_counts.get(node.name, 1)
        nstyle["size"] = min(size * scaling_factor, max_size)  # Scale and cap node size
        node.set_style(nstyle)

    # Improve tree spacing
    tree_style.branch_vertical_margin = 30  # Increase spacing for better visibility

    # Save the tree to a high-resolution PNG file
    output_file = "phylogenetic_tree.png"
    phylo_tree.render(output_file, w=3000, h=2000, tree_style=tree_style)

    print(f"Tree saved to {output_file}")

if __name__ == "__main__":
    main()
