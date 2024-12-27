import requests
import json

# Define the Pfam ID and base URL
pfam_id = "PF00151"
base_url = f"https://www.ebi.ac.uk/interpro/api/protein/reviewed/entry/pfam/{pfam_id}/"
headers = {"Accept": "application/json"}

# Initialize variables
results = []
page = 1

# Fetch the first page to get the total number of pages
response = requests.get(f"{base_url}?page=1", headers=headers)

if response.status_code == 200:
    total_pages = 5
    print(f"Total pages: {total_pages}")
else:
    print(f"Error: Unable to fetch data. Status code {response.status_code}")
    exit()

# Loop through all pages
while page <= total_pages:
    print(f"Fetching page {page}...")
    response = requests.get(f"{base_url}?page={page}", headers=headers)
    if response.status_code == 200:
        data = response.json().get("results", [])
        results.extend(data)
    else:
        print(f"Error fetching page {page}. Status code: {response.status_code}")
        break
    page += 1

# Save all results to a file
with open("all_results.json", "w") as f:
    json.dump(results, f, indent=4)

print("All results saved to all_results.json")
