import json

def extract_pfam_info(json_file):
    # Read the JSON file
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    # List to store the extracted information
    pfam_matches = []
    
    # Iterate through each result in the JSON data
    for result in data['results']:
        # Extract protein metadata
        protein_info = {
            'protein_name': result['metadata']['name'],
            'uniprot_id': result['metadata']['accession'],
            'organism': result['metadata']['source_organism']['scientificName']
        }
        
        # Extract PFAM domain information
        # We know there's only one entry because we queried for a specific PFAM domain
        pfam_entry = result['entries'][0]
        
        # Get the domain fragments (start and end positions)
        for location in pfam_entry['entry_protein_locations']:
            for fragment in location['fragments']:
                domain_info = {
                    **protein_info,  # Include all protein information
                    'domain_start': fragment['start'],
                    'domain_end': fragment['end'],
                    'domain_length': fragment['end'] - fragment['start'] + 1,
                    'protein_length': pfam_entry['protein_length'],
                    'score': location['score']
                }
                pfam_matches.append(domain_info)
    
    return pfam_matches

# Use the function to extract information
json_file = 'pfam_domain_positions.json'
matches = extract_pfam_info(json_file)

# Print the results in a formatted way
print("\nPFAM Domain Matches:")
print("-" * 80)
for match in matches:
    print(f"Protein: {match['protein_name']} ({match['uniprot_id']})")
    print(f"Organism: {match['organism']}")
    print(f"Domain position: {match['domain_start']}-{match['domain_end']} "
          f"(length: {match['domain_length']} aa)")
    print(f"Total protein length: {match['protein_length']} aa")
    print(f"Score: {match['score']}")
    print("-" * 80)

# Optional: Save to a more structured format like CSV for further analysis
import pandas as pd

# Convert the matches to a DataFrame
df = pd.DataFrame(matches)

# Save to CSV
df.to_csv('pfam_domain_positions.csv', index=False)
print("\nResults have been saved to 'pfam_domain_positions.csv'")