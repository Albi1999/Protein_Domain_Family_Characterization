# ONTOLOGY MARLON 

import requests
import pandas as pd
import xml.etree.ElementTree as ET
from scipy.stats import fisher_exact
from wordcloud import WordCloud
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
from statsmodels.stats.multitest import multipletests
from collections import defaultdict
from goatools import obo_parser




""" FUNCTION TASK """



# TODO : we basically did this above already for taxonomy task and just here neatly written into a function, so we could maybe just do it once in the whole code later on

# PRE-STEP 1 : Load Protein IDs of Family 
def load_protein_ids(psiblast_file, hmm_file, e_threshold=0.001):
    """Load protein IDs from PSI-BLAST and HMM search results."""
    psiblast_df = pd.read_csv(psiblast_file)
    hmm_df = pd.read_csv(hmm_file)
    
    filtered_hmm_proteins = hmm_df[hmm_df['E-value'] <= e_threshold]['uniprot_id']
    psiblast_proteins = set(psiblast_df['uniprot_id'])
    hmm_proteins = set(filtered_hmm_proteins)
    
    return list(psiblast_proteins.union(hmm_proteins))

# STEP 1 : For each Protein ID in our family, fetch its GO annotation 
def fetch_go_annotations(protein_id):
    """
    Fetch and categorize GO annotations for a given protein ID from the UniProt API.
    
    Args:
        protein_id (str): The UniProt ID of the protein
        
    Returns:
        dict: A dictionary containing:
            - Categorized GO terms separated by molecular function, biological process, 
              and cellular component (new format)
    """
    # Define the UniProt API URL for XML data
    url = f"https://rest.uniprot.org/uniprotkb/{protein_id}.xml"

    
    try:
        # Fetch the XML data from UniProt
        response = requests.get(url)
        response.raise_for_status()
        
        # Initialize our data structures
    
        categorized_terms = {
            'molecular_function': [],
            'biological_process': [],
            'cellular_component': []
        }
        
        # Set up namespace for XML parsing
        namespaces = {'ns': 'http://uniprot.org/uniprot'}
        root = ET.fromstring(response.content)
        
        # Find all GO term references in the XML
        for db_ref in root.findall(".//ns:dbReference[@type='GO']", namespaces):
            go_id = db_ref.attrib.get('id')
            term = db_ref.find("ns:property[@type='term']", namespaces)

            go_term = term.get('value')
    
     
            
            if go_id and term is not None:
                # Store in original format
                term_value = term.attrib['value']

                
                # Categorize based on prefix
                if term_value.startswith('F:'):
                    categorized_terms['molecular_function'].append({
                        'id': go_id,
                        'term': term_value[2:]  # Remove 'F:' prefix
                    })
                elif term_value.startswith('P:'):
                    categorized_terms['biological_process'].append({
                        'id': go_id,
                        'term': term_value[2:]  # Remove 'P:' prefix
                    })
                elif term_value.startswith('C:'):
                    categorized_terms['cellular_component'].append({
                        'id': go_id,
                        'term': term_value[2:]  # Remove 'C:' prefix
                    })
        
        return {
            'categorized': categorized_terms  # New categorized format
}
        
    except requests.exceptions.RequestException as e:
        print(f"Error fetching GO annotations for {protein_id}: {e}")
        return {
            'categorized': {
                'molecular_function': [],
                'biological_process': [],
                'cellular_component': []
            }
        }
    
# STEP 2 : GO ANNOTATIONS OF SWISSPROT .XML FILE 
# Let's add some debugging to help understand what's happening
# here we see that the big .xml file has the same structure as the small ones 
# we already analyzed ; thus,we can use the same parsing structure, but this time directly
# just collect the counts of GO terms, because that is all we need (no diff. categories, would just make our code slower)
def print_swissprot_file(swissprot_xml_path, length = 50):
    """
    Just to look at the first few lines to see the structure
    """

  
    with open(swissprot_xml_path, 'r') as f:
        print("First length lines of the file:")
        for i, line in enumerate(f):
            if i < length:
                print(line.strip())
            else:
                break
    



def parse_swissprot_go_terms(swissprot_xml_path, family_proteins):
    """
    Parse GO terms from SwissProt XML file, excluding proteins from our family.
    
    Args:
        swissprot_xml_path (str): Path to the SwissProt XML file
        family_proteins (set): Set of UniProt IDs in our protein family
    
    Returns:
        tuple: (go_term_counts dictionary, total proteins processed)
    """
    # Initialize counters
    go_term_counts = defaultdict(int)
    total_proteins = 0
    skipped_proteins = 0
    
    # Set up namespace for XML parsing
    namespaces = {'ns': 'http://uniprot.org/uniprot'}
    
    # Use iterparse for memory-efficient parsing
    context = ET.iterparse(swissprot_xml_path, events=('end',))
    
    print("Starting to parse SwissProt XML...")
    
    for event, elem in context:
        if elem.tag.endswith('entry'):
            # Get the UniProt ID for this protein
            accession = elem.find(".//ns:accession", namespaces)
            if accession is not None:
                uniprot_id = accession.text
                
                # Skip if this protein is in our family (we need this for the enrichment task to create the contigency table later on)
                if uniprot_id in family_proteins:
                    skipped_proteins += 1
                else:
                    # Process GO terms for non-family proteins
                    for db_ref in elem.findall(".//ns:dbReference[@type='GO']", namespaces):
                        go_id = db_ref.attrib.get('id') # the GO id 
                        if go_id:
                            go_term_counts[go_id] += 1
                    total_proteins += 1
            
            # Clear the element to save memory
            elem.clear()
            
            # Print progress periodically
            if (total_proteins + skipped_proteins) % 10000 == 0:
                print(f"Processed {total_proteins} proteins "
                      f"(skipped {skipped_proteins} family proteins)...")
              #  break
    
    return go_term_counts, total_proteins




def calculate_go_enrichment(my_go_counts, my_total_proteins, 
                          swissprot_go_counts, swissprot_total_proteins):
    results = []
    
    for go_id, my_count in my_go_counts.items():
        # Get count from SwissProt (i.e. how often was this GO id found in all of SwissProt with exception of
        # the Proteins found in the family?)
        swissprot_count = swissprot_go_counts.get(go_id, 0) # if isn't found, sets to count = 0 automatically 

        # Create the 2x2 contingency table for Fisher's exact test
        # The table looks like this:
        #                   Protein in family    Protein not in family (i.e. all in SwissProt - family proteins)
        # Has GO term            a                    b
        # No GO term             c                    d
        
        # Contingency table calculations:
        a = my_count  # Proteins with this GO term in family
        
    
        b = swissprot_count  # Proteins with GO term in SwissProt - the ones in Family
        
        c = my_total_proteins - a  # Proteins without GO term in family
   
        d = swissprot_total_proteins - b # Proteins that don't have the GO term in SwissProt - the ones in Family
        
        # Verify all values are non-negative before creating contingency table
        if all(x >= 0 for x in [a, b, c, d]): # TODO : remove 
            contingency_table = [[a, b], [c, d]]
            
            # Perform Fisher's exact test
            # We ask : is the GO term appearing more often in our family than we would expect by random chance ?
            # The null hypothesis (H0) is: "The proportion of proteins with this GO term in our family 
            # is the same as the proportion in the SwissProt dataset (without the protein in the family)." 
            # In other words, under H0, getting the GO term is independent of being in our family (so it doesn't represent the family)
            # Alternative Hypothesis (H1) depends on what tail to use 
            #Right-tail (greater): Our family has a higher proportion of this GO term than SwissProt
            #Left-tail (less): Our family has a lower proportion of this GO term than SwissProt
            #Two-tail (two-sided): The proportion is different (either higher or lower)
            #Fisher's exact test calculates the probability of seeing our observed data (or more extreme) under the null hypothesis.
            #A very small p-value (like < 0.05) tells us:
            #Two-tail: This GO term's frequency is significantly different from SwissProt
            #Right-tail: This GO term is significantly enriched in our family(overrepresented)
            #Left-tail: This GO term is significantly depleted in our family(underrepresented)

            odds_ratio, pvalue_two_tail = fisher_exact(contingency_table, alternative='two-sided')
            # TODO : including both the p-values for now, we have to understand when to use what (like asked in the task), 
            # TODO : i.e. how we ordered the confusion matrix (contingency table)
            _, pvalue_greater = fisher_exact(contingency_table, alternative='greater')
            _, pvalue_less = fisher_exact(contingency_table, alternative='less')
            
          
            my_proportion = my_count / my_total_proteins if my_total_proteins > 0 else 0
            swissprot_proportion = swissprot_count / swissprot_total_proteins if swissprot_total_proteins > 0 else 0
            # Fold Enrichment
            # TODO : see if the argumentation in the next comment makes sense (send email to prof)
            if swissprot_count == 0: # When the swissprot count is 0, it means that : 
                                     # When collecting the GO terms of SwissProt, we skipped over the proteins in our family
                                     # Thus, if no protein in SwissProt has this GO term, ONLY the protein in the family itself 
                                     # has that GO term (compared to ALL of SwissProt), thus in the WordCloud later on
                                     # we want to especially show the term of this GO id and will thus give it
                                     # 'inf' amount (infinite) for now
                if my_proportion > 0:
                    fold_enrichment = float('inf')
                else:
                    fold_enrichment = 0
            else:
                fold_enrichment = my_proportion/swissprot_proportion
       
     
            
            results.append({
                'GO_ID': go_id,
                'Count_Dataset': my_count,
                'Count_SwissProt': swissprot_count,
                'Percentage_Dataset': round(my_proportion * 100, 2),
                'Percentage_SwissProt': round(swissprot_proportion * 100, 10),
                'Fold_Enrichment': round(fold_enrichment,2),
                'P_Value_Two_Tail': pvalue_two_tail,
                'P_Value_Greater': pvalue_greater,
                'P_Value_Less': pvalue_less
            })
    

    
    # Convert to DataFrame and sort by p-value
    df_results = pd.DataFrame(results)
    if not df_results.empty:
        df_results = df_results.sort_values('P_Value_Two_Tail')

    df_results.to_csv("enrichment_results.csv")
    
    return df_results


# HELPER FUNCTION FOR STEP 2
def extract_go_terms_for_enrichment(protein_go_data):
    """
    Extract GO term counts from the protein annotation data.
    
    Args:
        protein_go_data (dict): Dictionary of protein annotations as provided
        
    Returns:
        tuple: (go_term_counts, total_proteins)
            - go_term_counts: Dictionary mapping GO terms to their counts
            - total_proteins: Total number of proteins in the dataset
    """
    go_term_counts = {}
    total_proteins = len(protein_go_data)
    
    # Iterate through each protein
    for protein_id, data in protein_go_data.items():
        # Get the categorized GO terms
        categories = data['categorized']
        
        # Process each category (molecular_function, biological_process, cellular_component)
        for category_terms in categories.values():
            # Count each GO term
            for term_info in category_terms:
                go_id = term_info['id']
                # Increment count if term exists, otherwise set to 1
                go_term_counts[go_id] = go_term_counts.get(go_id, 0) + 1
    
    return go_term_counts, total_proteins



# HELPER FUNCTION FOR STEP 3
# TODO : easier ? 
def create_go_id_to_term_mapping(family_data):
    """
    Creates a dictionary mapping GO IDs to their terms from the family data.
    Needed to create the word cloud based on the terms and not GO ID's
    
    Args:
        family_data (dict): Your dictionary containing protein annotations
        
    Returns:
        dict: Mapping of GO IDs to their terms
    """
    go_id_to_term = {}
    
    # Iterate through each protein's GO annotations
    for protein_id, annotations in family_data.items():
        for category in ['molecular_function', 'biological_process', 'cellular_component']:
            if 'categorized' in annotations and category in annotations['categorized']:
                for annotation in annotations['categorized'][category]:
                    go_id_to_term[annotation['id']] = annotation['term']
    
    return go_id_to_term








# STEP 4 
def analyze_go_hierarchy():
    # First, we downloaded the go.obo file so we can parse it 
    go_obo = obo_parser.GODag('go.obo')
    
    # Read our enrichment results
    df = pd.read_csv("enrichment_results.csv")
    
    # Filter for significantly enriched terms
    enriched_terms = df[
        (df['P_Value_Two_Tail'] < 0.05) &
        (df['P_Value_Greater'] < 0.05)
    ]
    
    # Create a dictionary to store branch information
    branch_info = {}
    
    # For each enriched term, traverse up its ancestry
    for _, row in enriched_terms.iterrows():
        go_id = row['GO_ID']
        if go_id in go_obo:
            term = go_obo[go_id]
            
            # Get all ancestors (parents) up to the root of the DAG (since we use get_all_parents we do that here! get_parents would just get the direct parents)
            ancestors = term.get_all_parents()
            
            # Add information about this term to all its ancestor branches
            for ancestor_id in ancestors:
                if ancestor_id not in branch_info:
                    branch_info[ancestor_id] = {
                        'term_name': go_obo[ancestor_id].name,
                        'enriched_children': [],
                        'total_significance': 0,
                        'depth': go_obo[ancestor_id].depth,
                    }

                # TODO : correct ????
                # Our go_id is a child to the current ancestors (note that this is not necessarily a direct child, but maybe also much more down in the tree somewhere)
                branch_info[ancestor_id]['enriched_children'].append({
                    'id': go_id,
                    'name': term.name,
                    'p_value': row['P_Value_Two_Tail']
                })
                # Add -log(p-value) to measure significance
                branch_info[ancestor_id]['total_significance'] += -np.log10(row['P_Value_Two_Tail'])
    
    # Filter for high-level terms (lower depth) with multiple enriched children
    significant_branches = {
        go_id: info for go_id, info in branch_info.items() # take each key,value of the branch_info dictionary
        if len(info['enriched_children']) >= 2  # At least 2 enriched children
        and info['depth'] <= 3  # High-level term (adjust this threshold as needed)
    }
    
    # Sort branches by their total significance
    sorted_branches = sorted(
        significant_branches.items(),
        key=lambda x: x[1]['total_significance'],
        reverse=True
    )
    
    # Create a list to store the branch information
    branch_data = []

    # Convert the branch information into a format suitable for a DataFrame
    for go_id, info in sorted_branches[:20]:  # Top 20 branches
        branch_data.append({
            'GO_ID': go_id,
            'Branch_Name': info['term_name'],
            'Hierarchy_Depth': info['depth'],
            'Number_Enriched_Terms': len(info['enriched_children']),
            'Total_Significance_Score': info['total_significance']
        })

    # Create a DataFrame and save to CSV
    branches_df = pd.DataFrame(branch_data)
    branches_df.to_csv('enriched_branches.csv', index=False)



def main():
    psiblast_file = "psiblast_parsed.csv"
    hmm_file = "hmmsearch_output.csv"
    protein_ids = load_protein_ids(psiblast_file, hmm_file)


    ######## STEP 1 ###########
    print("Fetching GO annotations...")
    family_annotations = {}
    for pid in tqdm(protein_ids, desc="Fetching GO annotations"):
        family_annotations[pid] = fetch_go_annotations(pid)


    print(family_annotations)

    go_counts_family, num_proteins_family = extract_go_terms_for_enrichment(family_annotations)

    print(go_counts_family)
    print(num_proteins_family)

       ######## STEP 1 END ###########

    # This step we already done and the files can be found here in the project (enrichment_results.csv)
    # Rerunning takes some time 
       ######## STEP 2 ###########
    """
    go_counts_swissprot, num_proteins_swissprot = parse_swissprot_go_terms("uniprot_sprot.xml", protein_ids)

    print(go_counts_swissprot)
    print(num_proteins_swissprot)
 

    _ = calculate_go_enrichment(go_counts_family, num_proteins_family,
                                            go_counts_swissprot, num_proteins_swissprot)
    """

       ######## STEP 2 END ###########

       ######## STEP 3 ###########
    # Read the enrichment results
    df = pd.read_csv("enrichment_results.csv")

    # Get the terms to the GO ids from the family data
    go_id_to_term = create_go_id_to_term_mapping(family_annotations)


    # Filter for significantly enriched terms based on two tail and right tail p-values
    enriched_terms = df[
    (df['P_Value_Two_Tail'] < 0.05) &
    (df['P_Value_Greater'] < 0.05)
    ]



    # Create word frequencies using the actual GO terms instead of IDs
    word_frequencies = {}
    for _, row in enriched_terms.iterrows():
        go_id = row['GO_ID']
        if go_id in go_id_to_term:  # Make sure we have the term for this ID
            term = go_id_to_term[go_id]
            # Use fold enrichment as weight, handling infinite values
            weight = 60000 if np.isinf(row['Fold_Enrichment']) else row['Fold_Enrichment'] 
            # TODO : we set to 60000 because we looked into the .csv and the highest Fold_Enrichment that 
            # did not have any

            word_frequencies[term] = weight

    # Create and display the word cloud
    wordcloud = WordCloud(
        width=1200, 
        height=800,
        background_color='white',
        prefer_horizontal=0.7,
        max_words=50,  # Limit to top 50 terms for better readability
        min_font_size=10,
        max_font_size=60
    ).generate_from_frequencies(word_frequencies)

    # Plot and save the word cloud
    plt.figure(figsize=(20, 12))
    plt.imshow(wordcloud, interpolation='bilinear')
    plt.axis('off')
    plt.title('GO Term Enrichment Word Cloud', fontsize=16, pad=20)
    plt.savefig('go_enrichment_wordcloud.png', dpi=300, bbox_inches='tight')
    plt.close()

    # Print out the enriched terms for verification
    print("\nTop enriched GO terms:")
    sorted_terms = sorted(word_frequencies.items(), key=lambda x: x[1], reverse=True)
    for term, weight in sorted_terms[:10]:
        print(f"\nTerm: {term}")
        print(f"Weight in word cloud: {weight:.2f}")

    
       ######## STEP 3 END ###########


    ######## STEP 4 ###########
    analyze_go_hierarchy()
    ######## STEP 4 END ###########


if __name__ == "__main__":
    main()



