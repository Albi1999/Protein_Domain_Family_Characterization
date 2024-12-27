import requests
import pandas as pd
import time
from tqdm import tqdm  # For progress bars

class InterProDomainFinder:
    def __init__(self):
        """Initialize the InterPro API client with base URL"""
        self.base_url = "https://www.ebi.ac.uk/interpro/api"
        
    def get_domain_positions(self, uniprot_id, pfam_id):
        """
        Get the positions of a specific Pfam domain in a protein
        
        Args:
            uniprot_id (str): UniProt accession (e.g., 'P12345')
            pfam_id (str): Pfam identifier (e.g., 'PF00151')
            
        Returns:
            list: List of dictionaries containing domain positions
        """
        # Construct the API URL for protein entries
        url = f"{self.base_url}/protein/reviewed/{uniprot_id}"
        
        # Make the API request
        response = requests.get(url, headers={'Accept': 'application/json'})
        
        if not response.ok:
            print(f"Error fetching data for {uniprot_id}: {response.status_code}")
            return []
            
        data = response.json()
        domain_positions = []
        
        # Extract entries containing domain information
        if 'entries' in data:
            for entry in data['entries']:
                # Look for Pfam entries matching our domain
                if entry.get('database') == 'pfam' and entry.get('accession') == pfam_id:
                    for location in entry.get('entry_protein_locations', []):
                        position = {
                            'uniprot_id': uniprot_id,
                            'start': location['fragments'][0]['start'],
                            'end': location['fragments'][0]['end'],
                            'score': location.get('score', None)
                        }
                        domain_positions.append(position)
                        
        return domain_positions

    def get_positions_for_protein_list(self, protein_list, pfam_id):
        """
        Get domain positions for a list of proteins
        
        Args:
            protein_list (list): List of UniProt accessions
            pfam_id (str): Pfam identifier
            
        Returns:
            DataFrame: Contains domain positions for all proteins
        """
        all_positions = []
        
        # Process each protein with a progress bar
        for protein in tqdm(protein_list, desc="Fetching domain positions"):
            # Add delay to respect API rate limits
            time.sleep(1)
            positions = self.get_domain_positions(protein, pfam_id)
            all_positions.extend(positions)
            
        return pd.DataFrame(all_positions)

# Example usage
if __name__ == "__main__":
    # Initialize the domain finder
    finder = InterProDomainFinder()
    
    # Example protein list and Pfam domain
    proteins = ["Q53H76", "Q8WWY8", "Q9Y5X9"]  # Your protein IDs
    pfam_id = "PF00151"  # Your Pfam domain
    
    # Get domain positions
    positions_df = finder.get_positions_for_protein_list(proteins, pfam_id)
    
    # Print results
    print("\nDomain positions found:")
    print(positions_df)
    
    # Save to CSV
    positions_df.to_csv('domain_positions.csv', index=False)
    print("\nResults saved to domain_positions.csv")