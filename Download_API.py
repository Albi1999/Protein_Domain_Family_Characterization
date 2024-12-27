import requests

# Define the Pfam ID and API URL
pfam_id = "PF00151"  # Replace with your Pfam ID
url = f"https://www.ebi.ac.uk/interpro/api/protein/reviewed/entry/pfam/{pfam_id}/"

# Perform the API request
response = requests.get(url, headers={"Accept": "application/json"})

if response.status_code == 200:
    # Save the response to a JSON file
    with open("results.json", "w") as f:
        f.write(response.text)
    print("Results saved to results.json")
else:
    print(f"Error: Unable to fetch data. Status code {response.status_code}")
