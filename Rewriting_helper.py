"""Used for rewriting the hmm_output and psiblast_output files such that they match the .csv structure of pfam_domain_positions.
This makes it easier to compare the .csv files directly (also without any code), so we can see if the results of metrics make sense
and how good our model is generally (based on e-scores and how much similar matches of proteins are found atleast)"""

# TODO : implement for hmm, but first we need to understand how we deal with it having multiple domain hits
import csv
import re 

def parse_psiblast_output(input_file):
    results = []
    
    with open(input_file, 'r') as f:
        for line in f:
            # Skip empty lines
            if not line.strip():
                continue
                
            # Split the line by tabs or multiple spaces
            parts = re.split(r'\s+', line.strip())
            
            if len(parts) >= 8:  # Make sure we have all required fields
                query_id = parts[0]
                subject_id = parts[1]
                
                # Extract UniProt ID and organism from subject_id
                # Format is usually sp|UniprotID|Name
                subject_parts = subject_id.split('|')
                if len(subject_parts) >= 2:
                    uniprot_id = subject_parts[1]
                    
                    # Create result dictionary
                    result = {
                        'protein_name': subject_id,
                        'uniprot_id': uniprot_id,
                        'organism': 'N/A',  # PSIBLAST output doesn't include organism
                        'domain_start': int(parts[4]),  # sstart
                        'domain_end': int(parts[5]),    # send
                        'domain_length': int(parts[5]) - int(parts[4]) + 1,
                        'E-value': float(parts[7])  
                    }
                    results.append(result)
    
    return results

def write_csv(results, output_file):
    if not results:
        return
    # Notice that we skip the start & end positions in the query domain (i.e. the PSSM here), as we are only interested in where we found matches in the sequence of SwissProt we looked through
    fieldnames = ['protein_name', 'uniprot_id', 'organism', 'domain_start', 
                 'domain_end', 'domain_length', 'E-value']
    
    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)


input_file = 'psiblast_search_output.txt'
output_file = 'psiblast_parsed.csv'

results = parse_psiblast_output(input_file)
write_csv(results, output_file)
print(f"Processed {len(results)} PSIBLAST matches")
