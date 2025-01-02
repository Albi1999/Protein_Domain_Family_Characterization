import requests
import pandas as pd
import xml.etree.ElementTree as ET
from scipy.stats import fisher_exact
from wordcloud import WordCloud
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
import random
import obonet
import networkx as nx
from statsmodels.stats.multitest import multipletests

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
# Step 2: Fetch GO Annotations from UniProt
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
# Step 3: Fetch Random Proteins
# -----------------------------------------------------------------------------

def fetch_random_proteins(batch_size=100, total_proteins=500):
    """Fetch a list of random reviewed UniProt protein IDs."""
    url = "https://rest.uniprot.org/uniprotkb/stream?query=reviewed:true&format=list"
    try:
        response = requests.get(url)
        response.raise_for_status()
        all_proteins = response.text.splitlines()
        selected_proteins = random.sample(all_proteins, min(total_proteins, len(all_proteins)))
        return [selected_proteins[i:i + batch_size] for i in range(0, len(selected_proteins), batch_size)]
    except requests.exceptions.RequestException as e:
        print(f"Error fetching random proteins: {e}")
        return []

# -----------------------------------------------------------------------------
# Step 4: Load GO Ontology
# -----------------------------------------------------------------------------

def load_go_ontology():
    """Load the GO ontology from OBO file."""
    url = "http://purl.obolibrary.org/obo/go/go-basic.obo"
    graph = obonet.read_obo(url)
    return graph

# -----------------------------------------------------------------------------
# Step 5: Flatten Annotations for Analysis
# -----------------------------------------------------------------------------

def flatten_annotations(annotation_dict):
    """Flatten GO annotations into a list of GO terms."""
    flat_terms = []
    for annotations in annotation_dict.values():
        flat_terms.extend([a["GO_ID"] for a in annotations])
    return flat_terms

# -----------------------------------------------------------------------------
# Step 6: Enrichment Analysis
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
# Step 7: Visualize Enriched Terms
# -----------------------------------------------------------------------------

def plot_wordcloud(enrichment_results, annotations):
    """Generate and save a word cloud of enriched GO terms."""
    term_names = {a["GO_ID"]: a["Term"] for ann_list in annotations.values() for a in ann_list}
    enriched_with_names = {term_names[go_id]: -np.log10(p) for go_id, p in enrichment_results.items() if go_id in term_names}

    if not enriched_with_names:
        print("No enriched terms found. Word cloud will not be generated.")
        return

    wordcloud = WordCloud(width=1000, height=600, background_color="white", colormap="viridis").generate_from_frequencies(enriched_with_names)
    wordcloud.to_file("enriched_terms_wordcloud.png")
    print("Word cloud saved as enriched_terms_wordcloud.png")

    plt.figure(figsize=(12, 6))
    plt.imshow(wordcloud, interpolation="bilinear")
    plt.axis("off")
    plt.title("Word Cloud of Enriched GO Terms")
    plt.tight_layout()
    plt.show()

def plot_branch_enrichment(enrichment_results, go_graph):
    """Plot GO branch enrichment."""
    branch_scores = {}
    for go_id, p_value in enrichment_results.items():
        try:
            parents = nx.ancestors(go_graph, go_id)
            for parent in parents:
                if parent not in branch_scores:
                    branch_scores[parent] = []
                branch_scores[parent].append(p_value)
        except nx.NetworkXError:
            continue

    significant_branches = {branch: np.mean(scores)
                            for branch, scores in branch_scores.items() if len(scores) >= 3}

    # Limit to top 50 branches by score
    sorted_branches = sorted(significant_branches.items(), key=lambda x: x[1], reverse=True)[:50]
    branches, scores = zip(*sorted_branches)

    plt.figure(figsize=(14, 10))
    plt.barh(range(len(branches)), scores, color="steelblue")
    plt.yticks(range(len(branches)), [go_graph.nodes[branch]['name'] for branch in branches], fontsize=8)
    plt.xlabel('Mean p-value')
    plt.title('Top 50 GO Branch Enrichments')
    plt.tight_layout()
    plt.savefig("go_enrichment_branches.png", dpi=300, bbox_inches="tight")
    plt.legend(["Top 50 Branch Enrichments"], loc="lower right")
    print("GO branch enrichment plot saved as go_enrichment_branches.png")

# -----------------------------------------------------------------------------
# Step 8: Write Summary and Results to File
# -----------------------------------------------------------------------------

def write_summary_and_results(enrichment_results, family_annotations, go_graph):
    """Write a summary and detailed results to a text file."""
    with open("enrichment_results.txt", "w") as f:
        # Write summary
        f.write("SUMMARY\n")
        f.write("========\n")
        f.write(f"Number of enriched GO terms: {len(enrichment_results)}\n")
        f.write(f"Top enriched term: {max(enrichment_results, key=enrichment_results.get, default='None')}\n")
        f.write("\n\n")

        # Write detailed results
        f.write("DETAILED RESULTS\n")
        f.write("================\n")
        for go_id, p_value in enrichment_results.items():
            term_name = next((a["Term"] for ann_list in family_annotations.values() 
                             for a in ann_list if a["GO_ID"] == go_id), go_id)
            f.write(f"{go_id}: {term_name} (p-value: {p_value:.2e})\n")

        # Write branch scores
        f.write("\n\nSIGNIFICANT GO BRANCHES\n")
        f.write("========================\n")
        branch_scores = {}
        for go_id, p_value in enrichment_results.items():
            try:
                parents = nx.ancestors(go_graph, go_id)
                for parent in parents:
                    if parent not in branch_scores:
                        branch_scores[parent] = []
                    branch_scores[parent].append(p_value)
            except nx.NetworkXError:
                continue

        significant_branches = {branch: np.mean(scores) 
                                for branch, scores in branch_scores.items() if len(scores) >= 3}
        for branch, score in sorted(significant_branches.items(), key=lambda x: x[1]):
            branch_name = go_graph.nodes[branch].get('name', branch)
            f.write(f"{branch}: {branch_name} (mean p-value: {score:.2e})\n")

# -----------------------------------------------------------------------------
# Main Script
# -----------------------------------------------------------------------------

def main():
    print("Loading GO ontology...")
    go_graph = load_go_ontology()

    psiblast_file = "psiblast_parsed.csv"
    hmm_file = "hmmsearch_output.csv"
    protein_ids = load_protein_ids(psiblast_file, hmm_file)

    print("Fetching GO annotations...")
    family_annotations = {}
    for pid in tqdm(protein_ids, desc="Fetching GO annotations"):
        family_annotations[pid] = fetch_go_annotations(pid)

    print("Fetching background annotations...")
    background_annotations = {}
    background_batches = fetch_random_proteins(batch_size=50, total_proteins=500)
    for batch in tqdm(background_batches, desc="Processing background proteins"):
        for pid in batch:
            background_annotations[pid] = fetch_go_annotations(pid)

    print("Calculating enrichment...")
    family_terms = flatten_annotations(family_annotations)
    background_terms = flatten_annotations(background_annotations)

    unique_go_terms = set(family_terms)
    enrichment_results = {}
    pvalues = []
    terms = []

    for term in unique_go_terms:
        _, p_value = fisher_exact([
            [family_terms.count(term), len(family_terms) - family_terms.count(term)],
            [background_terms.count(term), len(background_terms) - background_terms.count(term)]
        ], alternative='greater')

        pvalues.append(p_value)
        terms.append(term)

    rejected, p_corrected, _, _ = multipletests(pvalues, method='fdr_bh')

    for term, p_value, significant in zip(terms, p_corrected, rejected):
        if significant:
            enrichment_results[term] = p_value

    print("Generating visualizations...")
    plot_wordcloud(enrichment_results, family_annotations)
    plot_branch_enrichment(enrichment_results, go_graph)

    print("Writing results to file...")
    write_summary_and_results(enrichment_results, family_annotations, go_graph)

if __name__ == "__main__":
    main()
