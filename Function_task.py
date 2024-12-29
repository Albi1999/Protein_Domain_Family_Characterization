# Function Task: Perform Gene Ontology (GO) enrichment analysis for a set of proteins.

import requests
import xml.etree.ElementTree as ET
from scipy.stats import fisher_exact
from wordcloud import WordCloud
import matplotlib.pyplot as plt
import random
import numpy as np

# -----------------------------------------------------------------------------
# Step 1: Fetch GO Annotations from UniProt
# -----------------------------------------------------------------------------

def fetch_go_annotations(protein_id):
    """Fetch GO annotations for a given protein ID from the UniProt API."""
    url = f"https://rest.uniprot.org/uniprotkb/{protein_id}.xml"
    try:
        response = requests.get(url)
        response.raise_for_status()
        go_terms = []
        namespaces = {'ns': 'http://uniprot.org/uniprot'}
        root = ET.fromstring(response.content)
        for db_ref in root.findall(".//ns:dbReference[@type='GO']", namespaces):
            go_id = db_ref.attrib.get('id')
            term = db_ref.find("ns:property[@type='term']", namespaces)
            if go_id and term is not None:
                go_terms.append({"GO_ID": go_id, "Term": term.attrib['value']})
        return go_terms
    except requests.exceptions.RequestException as e:
        print(f"Error fetching GO annotations for {protein_id}: {e}")
        return []

# -----------------------------------------------------------------------------
# Step 2: Read Protein IDs from a File
# -----------------------------------------------------------------------------

def read_protein_ids(file_path):
    """Read protein IDs from a text file."""
    with open(file_path, "r") as f:
        return [line.strip() for line in f.readlines()]

# -----------------------------------------------------------------------------
# Step 3: Fetch Background GO Annotations
# -----------------------------------------------------------------------------

def fetch_random_proteins(batch_size=100, total_proteins=500):
    """Fetch a list of random reviewed UniProt protein IDs."""
    url = "https://rest.uniprot.org/uniprotkb/stream?query=reviewed:true&format=list"
    try:
        response = requests.get(url)
        response.raise_for_status()
        all_proteins = response.text.splitlines()
        selected_proteins = random.sample(all_proteins, total_proteins)
        return [selected_proteins[i:i + batch_size] for i in range(0, len(selected_proteins), batch_size)]
    except requests.exceptions.RequestException as e:
        print(f"Error fetching random proteins: {e}")
        return []

def fetch_background_annotations(protein_batches):
    """Fetch GO annotations for batches of proteins."""
    background_annotations = {}
    for batch in protein_batches:
        for pid in batch:
            background_annotations[pid] = fetch_go_annotations(pid)
    return background_annotations

# -----------------------------------------------------------------------------
# Step 4: Flatten Annotations for Analysis
# -----------------------------------------------------------------------------

def flatten_annotations(annotation_dict):
    """Flatten GO annotations into a list of GO terms."""
    flat_terms = []
    for annotations in annotation_dict.values():
        flat_terms.extend([a["GO_ID"] for a in annotations])
    return flat_terms

# -----------------------------------------------------------------------------
# Step 5: Enrichment Analysis
# -----------------------------------------------------------------------------

def calculate_enrichment(go_term, family_terms, background_terms):
    """Calculate enrichment of a GO term using Fisher's exact test."""
    family_count = family_terms.count(go_term)
    family_not = len(family_terms) - family_count
    background_count = background_terms.count(go_term)
    background_not = len(background_terms) - background_count

    contingency_table = [[family_count, background_count],
                         [family_not, background_not]]
    _, p_value = fisher_exact(contingency_table, alternative='greater')
    return p_value

# -----------------------------------------------------------------------------
# Step 6: Visualize Enriched Terms
# -----------------------------------------------------------------------------

def plot_wordcloud(enrichment_results, annotations, output_file="enriched_terms_wordcloud.png"):
    """Generate and save a word cloud of enriched GO terms."""
    term_names = {a["GO_ID"]: a["Term"] for ann_list in annotations.values() for a in ann_list}
    enriched_with_names = {term_names[go_id]: -np.log10(p) for go_id, p in enrichment_results.items() if go_id in term_names}

    if not enriched_with_names:
        print("No enriched terms found. Word cloud will not be generated.")
        return

    wordcloud = WordCloud(width=800, height=400, background_color="white").generate_from_frequencies(enriched_with_names)
    wordcloud.to_file(output_file)
    print(f"Word cloud saved as {output_file}")

    plt.figure(figsize=(10, 5))
    plt.imshow(wordcloud, interpolation="bilinear")
    plt.axis("off")
    plt.title("Enriched GO Terms")
    plt.tight_layout()
    plt.show()

# -----------------------------------------------------------------------------
# Main Script
# -----------------------------------------------------------------------------

def main():
    # Load family protein IDs
    protein_ids = read_protein_ids("cleaned_protein_ids.txt")
    family_annotations = {pid: fetch_go_annotations(pid) for pid in protein_ids}

    # Fetch background annotations in batches
    protein_batches = fetch_random_proteins(batch_size=50, total_proteins=200)
    background_annotations = fetch_background_annotations(protein_batches)

    # Flatten annotations for Fisher's test
    family_terms = flatten_annotations(family_annotations)
    background_terms = flatten_annotations(background_annotations)

    # Perform enrichment analysis
    unique_go_terms = set(family_terms)
    enrichment_results = {}
    for term in unique_go_terms:
        p_value = calculate_enrichment(term, family_terms, background_terms)
        if p_value < 0.05:  # Adjust threshold if needed
            enrichment_results[term] = p_value

    # Generate and save the word cloud
    plot_wordcloud(enrichment_results, family_annotations)

if __name__ == "__main__":
    main()
