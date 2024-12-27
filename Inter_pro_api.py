import requests

# Define the Pfam ID and API URL
pfam_id = "PF00151"  # Replace with your Pfam ID
url = f"https://www.ebi.ac.uk/interpro/api/protein/reviewed/entry/pfam/{pfam_id}/"


# Perform the API request
response = requests.get(url, headers={"Accept": "application/json"})

if response.status_code == 200:
    # Parse the JSON response
    proteins = response.json().get("results", [])
    for protein in proteins:
        protein_id = protein.get("metadata", {}).get("accession", "Unknown")
        locations = protein.get("entry_protein_locations", [])
        if locations:
            for loc in locations:
                fragments = loc.get("fragments", [])
                for fragment in fragments:
                    start = fragment.get("start", "Unknown")
                    end = fragment.get("end", "Unknown")
                    print(f"Protein: {protein_id}, Start: {start}, End: {end}")
        else:
            print(f"Protein: {protein_id}, Domain positions not available")
else:
    print(f"Error: Unable to fetch data. Status code {response.status_code}")
