# Extract UniProt IDs from a FASTA file
fasta_file = "trimmed_alignment.fasta"  # Update the path if needed
cleaned_protein_ids = []

with open(fasta_file, "r") as file:
    for line in file:
        if line.startswith(">"):  # Header line in FASTA
            # Extract the first token after '>'
            protein_id = line.split()[0][1:]  # Remove '>'
            # Keep only the second field (after the first '|')
            if "|" in protein_id:
                protein_id = protein_id.split("|")[1]
            cleaned_protein_ids.append(protein_id)

print(f"Cleaned Protein IDs: {cleaned_protein_ids}")


# Save the cleaned IDs to a text file
with open("cleaned_protein_ids.txt", "w") as output_file:
    output_file.write("\n".join(cleaned_protein_ids))
