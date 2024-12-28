import requests

def fetch_family_sequences_from_file(input_file, output_file="family_sequences.xml"):
    with open(input_file, "r") as file:
        protein_ids = [line.strip() for line in file.readlines()]
    
    base_url = "https://rest.uniprot.org/uniprotkb/stream"
    params = {
        "format": "xml",
        "query": " OR ".join([f"accession:{pid}" for pid in protein_ids])
    }
    response = requests.get(base_url, params=params)
    
    if response.status_code == 200:
        with open(output_file, "wb") as file:
            file.write(response.content)
        print(f"File saved as {output_file}")
    else:
        print(f"Failed to fetch data: {response.status_code}, {response.text}")

# Example usage
input_file = "cleaned_protein_ids.txt"  # Replace with the name of your txt file
fetch_family_sequences_from_file(input_file)
