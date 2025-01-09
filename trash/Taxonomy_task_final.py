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

# Process taxonomy data
def process_taxonomy(data, correct_column_name):
    taxonomy_dict = {}
    frequency_counts = {}
    
    for _, row in data.iterrows():
        lineage = row[correct_column_name].split(" > ")
        current = taxonomy_dict
        # Track the full path to maintain hierarchy information
        current_path = [] # such that we count occurences of terms in the correct "level" where they appear (i.e. always count just in the "column" of the linage)
        
        for level in lineage:
            current_path.append(level)
            path_key = " > ".join(current_path)
            
            # Count frequencies using the full path as key
            if path_key not in frequency_counts:
                frequency_counts[path_key] = 0
            frequency_counts[path_key] += 1
            
            if level not in current:
                current[level] = {}
            current = current[level]
    
    return taxonomy_dict, frequency_counts

# Create a Newick string for the taxonomy tree
def dict_to_newick(d, parent_abundance=None):
    newick = ""
    for key, sub_dict in d.items():
        size = parent_abundance.get(key, 1) if parent_abundance else 1
        sub_tree = dict_to_newick(sub_dict, parent_abundance)
        newick += f"({sub_tree}){key}:{size}," if sub_tree else f"{key}:{size},"
    return newick.rstrip(",")



# Fetch taxonomy data
def main():
    psiblast_file = "psiblast_parsed.csv"
    hmm_file = "hmmsearch_output.csv"
    protein_ids = load_protein_ids(psiblast_file, hmm_file)

    analyzer = TaxonomyAnalyzer()
    taxonomy_data = analyzer.fetch_taxonomy_info(protein_ids, "taxonomy_info.csv")

    print("Taxonomy file saved to: taxonomy_info.csv")

    # Correct column name
    correct_column_name = "Lineage"  # Use the correct column name

    # Create a nested dictionary of taxonomy
    taxonomy_dict, frequency_counts = process_taxonomy(taxonomy_data, correct_column_name)

    # Count relative abundance (of the different paths ! ; right now we don't really use that)
    abundance_counts = taxonomy_data[correct_column_name].value_counts().to_dict()

    

# TODO : abundance counts used here, but it doesn't show at all in the graph ; we need to ask professor if what we have now is already enough, then we should remove 
# TODO : this part with abundance counts
    newick_tree = f"({dict_to_newick(taxonomy_dict, abundance_counts)});"

    # Plot using ETE Toolkit
    phylo_tree = Tree(newick_tree, format=1)
    tree_style = TreeStyle()
    tree_style.show_leaf_name = False


    # Adjust node sizes (normalize and refine scaling)
    max_size = 50  # Increase max size for better differentiation
    scaling_factor = 2  # Further refine scaling for visual contrast
    for node in phylo_tree.traverse():
        # Get the full path from root to this node
        path = []
        current = node
        while current:
            if current.name:  # Skip empty names
                path.insert(0, current.name)
            current = current.up
        
        path_key = " > ".join(path)
        count = frequency_counts.get(path_key, 1)
        nstyle = NodeStyle()
        size = abundance_counts.get(node.name, 1)
        nstyle["size"] = min(size * scaling_factor, max_size)  # Scale and cap node size
        node.set_style(nstyle)
        # Add label with name and count
        node.add_face(TextFace(f"{node.name} ({count})", fsize=10), column=0)

    # Improve tree spacing
    tree_style.branch_vertical_margin = 30  # Increase spacing for better visibility

    # Save the tree to a high-resolution PNG file
    output_file = "phylogenetic_tree_freq.png"
    phylo_tree.render(output_file, w=3000, h=2000, tree_style=tree_style)

    print(f"Tree saved to: {output_file}")

if __name__ == "__main__":
    main()
