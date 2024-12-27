import json
import requests
import time
from typing import Dict, List, Optional
from collections import defaultdict

class InterProAPIFetcher:
    def __init__(self, base_url: str):
        self.base_url = base_url
        self.processed_count = 0
        self.all_results = []
        self.seen_accessions = set()  # Track unique protein accessions
        self.duplicate_count = 0

    def fetch_page(self, url: str) -> Optional[Dict]:
        max_retries = 3
        retry_delay = 2
        
        for attempt in range(max_retries):
            try:
                response = requests.get(url)
                response.raise_for_status()
                return response.json()
            except requests.exceptions.RequestException as e:
                print(f"Attempt {attempt + 1} failed: {str(e)}")
                if attempt < max_retries - 1:
                    print(f"Waiting {retry_delay} seconds before retrying...")
                    time.sleep(retry_delay)
                    retry_delay *= 2
                else:
                    print("Max retries reached. Moving on...")
                    return None

    def fetch_all_pages(self) -> List[Dict]:
        next_url = self.base_url
        total_count = None
        page_number = 1
        
        while next_url:
            print(f"\nFetching page {page_number}...")
            print(f"URL: {next_url}")
            
            page_data = self.fetch_page(next_url)
            
            if page_data is None:
                print("Failed to fetch page. Stopping pagination.")
                break
            
            if total_count is None:
                total_count = page_data['count']
                print(f"API reports total count: {total_count}")
            
            # Check for duplicates in this page
            new_proteins = []
            page_duplicates = 0
            
            for protein in page_data['results']:
                accession = protein['metadata']['accession']
                if accession in self.seen_accessions:
                    page_duplicates += 1
                    self.duplicate_count += 1
                else:
                    self.seen_accessions.add(accession)
                    new_proteins.append(protein)
            
            print(f"Page {page_number} stats:")
            print(f"- Proteins in response: {len(page_data['results'])}")
            print(f"- New unique proteins: {len(new_proteins)}")
            print(f"- Duplicates found: {page_duplicates}")
            
            self.all_results.extend(new_proteins)
            self.processed_count = len(self.all_results)
            
            next_url = page_data.get('next')
            page_number += 1
            
            time.sleep(1)
        
        print(f"\nFinal Statistics:")
        print(f"Total unique proteins: {len(self.all_results)}")
        print(f"Total duplicates found: {self.duplicate_count}")
        print(f"Total processed entries: {self.processed_count + self.duplicate_count}")
        
        return self.all_results

    def save_results(self, filename: str):
        output_data = {
            'count': len(self.all_results),
            'results': self.all_results
        }
        
        with open(filename, 'w') as f:
            json.dump(output_data, f, indent=2)
        print(f"\nResults saved to {filename}")
        print(f"File contains {len(self.all_results)} unique proteins")

def main():
    base_url = "https://www.ebi.ac.uk/interpro/api/protein/reviewed/entry/pfam/PF00151/"
    
    fetcher = InterProAPIFetcher(base_url)
    fetcher.fetch_all_pages()
    fetcher.save_results('complete_pfam_matches.json')

if __name__ == "__main__":
    main()