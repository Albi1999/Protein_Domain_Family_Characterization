import requests
import pandas as pd
from scipy.stats import fisher_exact
from wordcloud import WordCloud
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
import obonet
from statsmodels.stats.multitest import multipletests
import xml.etree.ElementTree as ET
from io import StringIO
import networkx as nx

# -----------------------------------------------------------------------------
# Step 1: Load Protein IDs
# -----------------------------------------------------------------------------
def load_protein_ids(psiblast_file, hmm_file, e_threshold=0.001):
    """Load protein IDs from PSI-BLAST and HMM search results."""
    psiblast_df = pd.read_csv(psiblast_file)
    hmm_df = pd.read_csv(hmm_file)
    
    filtered_hmm_proteins = hmm_df[hmm_df['E-value'] <= e_threshold]['uniprot_id']
    psiblast_proteins = set(psiblast_df['uniprot_id'])
    hmm_proteins = set(filtered_hmm_proteins)
    
    return list(psiblast_proteins.union(hmm_proteins))

# -----------------------------------------------------------------------------
# Step 2: Fetch GO Annotations
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
            category = db_ref.find("ns:property[@type='category']", namespaces)
            if go_id and term is not None:
                go_terms.append({
                    "GO_ID": go_id,
                    "Term": term.attrib['value'],
                    "Category": category.attrib['value'] if category is not None else "Unknown"
                })
        return go_terms
    except requests.exceptions.RequestException as e:
        print(f"Error fetching GO annotations for {protein_id}: {e}")
        return []

# -----------------------------------------------------------------------------
# Step 3: Fetch Background GO Annotations
# -----------------------------------------------------------------------------
def fetch_background_annotations(batch_size=100):
    """Fetch a list of random UniProt reviewed protein GO annotations."""
    url = "https://rest.uniprot.org/uniprotkb/stream?query=reviewed:true&format=tsv&fields=accession,id,go" 
    try:
        response = requests.get(url)
        response.raise_for_status()
        background_annotations = {}

        # Parse the response for GO annotations
        lines = response.text.splitlines()[1:]  # Skip header
        total_lines = len(lines)
        with tqdm(total=total_lines, desc="Parsing background annotations", unit="line") as pbar:
            for line in lines:
                fields = line.split("\t")
                protein_id = fields[0]
                go_annotations = []

                if len(fields) > 2 and fields[2]:
                    go_annotations = [term.split(";")[0] for term in fields[2].split("|")]

                background_annotations[protein_id] = go_annotations
                pbar.update(1)

        return background_annotations

    except requests.exceptions.RequestException as e:
        print(f"Error fetching background GO annotations: {e}")
        return {}

# -----------------------------------------------------------------------------
# Step 4: Enrichment Analysis
# -----------------------------------------------------------------------------
def calculate_enrichment(go_term, family_terms, background_terms):
    """Calculate enrichment of a GO term using Fisher's exact test."""
    family_count = family_terms.count(go_term)
    family_not = len(family_terms) - family_count
    background_count = background_terms.count(go_term)
    background_not = len(background_terms) - background_count

    contingency_table = [[family_count, background_count],
                         [family_not, background_not]]
    two_tail, right_tail = fisher_exact(contingency_table, alternative='two-sided'), fisher_exact(contingency_table, alternative='greater')
    return two_tail[1], right_tail[1]

# -----------------------------------------------------------------------------
# Step 5: Flatten Annotations for Analysis
# -----------------------------------------------------------------------------
def flatten_annotations(annotation_dict):
    """Flatten GO annotations into a list of GO IDs only."""
    flat_terms = []
    for key, annotations in annotation_dict.items():
        if isinstance(annotations, list) and annotations:  # Check if the list is non-empty
            if isinstance(annotations[0], dict):  # Family annotations: list of dicts with "GO_ID"
                flat_terms.extend([annotation["GO_ID"] for annotation in annotations])
            else:  # Background annotations: list of GO IDs
                flat_terms.extend(annotations)
    return flat_terms

# -----------------------------------------------------------------------------
# Step 6: Fetch GO Ontology
# -----------------------------------------------------------------------------
def fetch_go_ontology():
    """Fetch the GO ontology dynamically using the Gene Ontology PURL."""
    url = "http://purl.obolibrary.org/obo/go/go-basic.obo"
    try:
        response = requests.get(url)
        response.raise_for_status()
        go_graph = obonet.read_obo(StringIO(response.text))
        return go_graph
    except requests.exceptions.RequestException as e:
        print(f"Error fetching GO ontology: {e}")
        return None

# -----------------------------------------------------------------------------
# Step 7: Visualization - Word Cloud and Branch Enrichment
# -----------------------------------------------------------------------------
def plot_wordcloud(enrichment_results, go_annotations):
    """Generate and save a word cloud of enriched GO terms."""
    # Map GO IDs to their term names
    term_mapping = {term["GO_ID"]: term["Term"] for terms in go_annotations.values() for term in terms}

    # Replace GO IDs with term names and calculate scores
    term_scores = {
        term_mapping.get(go_id, go_id): -np.log10(p)
        for go_id, p in enrichment_results.items()
        if go_id in term_mapping
    }

    # Limit to the top 50 terms for better readability
    top_terms = dict(sorted(term_scores.items(), key=lambda item: item[1], reverse=True)[:50])

    # Create the word cloud with enhanced settings
    wordcloud = WordCloud(
        width=1200,
        height=800,
        background_color='white',
        prefer_horizontal=0.95,  # Increased horizontal preference
        max_words=40,
        min_font_size=10,
        max_font_size=60,
        colormap='viridis'
    ).generate_from_frequencies(top_terms)

    # Save and display the word cloud
    wordcloud.to_file("enriched_terms_wordcloud_final.png")
    plt.figure(figsize=(15, 10))
    plt.imshow(wordcloud, interpolation="bilinear")
    plt.axis("off")
    plt.title("GO Term Enrichment Word Cloud", fontsize=20)
    plt.tight_layout()
    plt.show()


# Branch enrichment plotting remains unchanged
def plot_branch_enrichment(enrichment_results, go_graph):
    """Generate a plot for branch enrichment."""
    branch_scores = {}
    for go_id, p_value in enrichment_results.items():
        if go_id not in go_graph.nodes:
            continue

        parents = nx.ancestors(go_graph, go_id)
        for parent in parents:
            if parent not in branch_scores:
                branch_scores[parent] = []
            branch_scores[parent].append(p_value)

    significant_branches = {branch: np.mean(scores) for branch, scores in branch_scores.items() if len(scores) >= 3}
    sorted_branches = sorted(significant_branches.items(), key=lambda x: x[1])[:50]

    branch_ids, scores = zip(*sorted_branches)
    branch_names = [go_graph.nodes[branch].get("name", branch) for branch in branch_ids]

    plt.figure(figsize=(14, 10))
    plt.barh(branch_names, scores, color="steelblue")
    plt.xlabel("Mean p-value")
    plt.ylabel("GO Branch")
    plt.title("Top Enriched GO Branches")
    plt.tight_layout()
    plt.savefig("go_enrichment_branches_final.png")
    plt.close()

# -----------------------------------------------------------------------------
# Step 8: Save Results
# -----------------------------------------------------------------------------
def save_enrichment_results(enrichment_results, go_annotations, family_terms, background_terms):
    """Save enrichment results to a CSV file with detailed information."""
    term_mapping = {term["GO_ID"]: term for terms in go_annotations.values() for term in terms}

    detailed_results = []
    for go_id, p_value in enrichment_results.items():
        # Fetch term details
        term_details = term_mapping.get(go_id, {})
        term_name = term_details.get("Term", "Unknown")
        category = term_details.get("Category", "Unknown")

        # Count occurrences in family and background
        family_count = family_terms.count(go_id)
        background_count = background_terms.count(go_id)

        # Add detailed data
        detailed_results.append({
            "GO_ID": go_id,
            "Term_Name": term_name,
            "Category": category,
            "Family_Count": family_count,
            "Background_Count": background_count,
            "Two_Tail_P": p_value,
            "Right_Tail_P": p_value  # If both two-tailed and right-tailed are the same
        })

    # Convert to DataFrame and save to CSV
    df = pd.DataFrame(detailed_results)
    df.to_csv("enrichment_results_final.csv", index=False)
    print("Enrichment results saved as enrichment_results_final.csv")

def save_summary(enrichment_results, family_terms, background_terms):
    """Save a detailed summary of the enrichment analysis."""
    with open("enrichment_summary_final.txt", "w") as f:
        f.write("Enrichment Analysis Summary\n")
        f.write("===========================\n")
        f.write(f"Total enriched terms: {len(enrichment_results)}\n")
        f.write(f"Total family GO terms: {len(set(family_terms))}\n")
        f.write(f"Total background GO terms: {len(set(background_terms))}\n")
        if enrichment_results:
            top_term = min(enrichment_results, key=enrichment_results.get)
            f.write(f"Most significant term: {top_term} with p-value {enrichment_results[top_term]:.2e}\n")
        f.write("\nAdditional Notes:\n")
        f.write("- Terms skipped due to missing ontology nodes are automatically excluded from visualizations.\n")

# -----------------------------------------------------------------------------
# Main Script
# -----------------------------------------------------------------------------
def main():
    psiblast_file = "psiblast_parsed.csv"
    hmm_file = "hmmsearch_output.csv"

    # Load family protein IDs
    protein_ids = load_protein_ids(psiblast_file, hmm_file)

    # Fetch GO annotations for family
    family_annotations = {}
    with tqdm(total=len(protein_ids), desc="Fetching GO annotations") as pbar:
        for pid in protein_ids:
            family_annotations[pid] = fetch_go_annotations(pid)
            pbar.update(1)

    # Fetch background GO annotations
    background_annotations = fetch_background_annotations()

    # Flatten annotations
    family_terms = flatten_annotations(family_annotations)
    background_terms = flatten_annotations(background_annotations)

    # Fetch GO ontology
    go_graph = fetch_go_ontology()
    if go_graph is None:
        print("Failed to fetch GO ontology. Exiting...")
        return

    # Enrichment analysis
    unique_go_terms = set(family_terms)
    enrichment_results = {}
    for term in tqdm(unique_go_terms, desc="Analyzing GO terms"):
        two_tail_p, right_tail_p = calculate_enrichment(term, family_terms, background_terms)
        if right_tail_p < 0.05 or two_tail_p < 0.05:
            enrichment_results[term] = two_tail_p

    # Save results and visualizations
    save_enrichment_results(enrichment_results, family_annotations, family_terms, background_terms)
    save_summary(enrichment_results, family_terms, background_terms)
    plot_wordcloud(enrichment_results, family_annotations)
    plot_branch_enrichment(enrichment_results, go_graph)

if __name__ == "__main__":
    main()
