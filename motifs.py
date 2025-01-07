import re
import csv
from Bio import SeqIO
import ast
import pandas as pd

def load_protein_ids(psiblast_file, hmm_file, e_threshold=0.001):
    """Load protein IDs from PSI-BLAST and HMM search results."""
    psiblast_df = pd.read_csv(psiblast_file)
    hmm_df = pd.read_csv(hmm_file)
    
    filtered_hmm_proteins = hmm_df[hmm_df['E-value'] <= e_threshold]['uniprot_id']
    print(filtered_hmm_proteins)
    psiblast_proteins = set(psiblast_df['uniprot_id'])
    hmm_proteins = set(filtered_hmm_proteins)
    
    return psiblast_proteins.union(hmm_proteins)

def load_disordered_regions(mobidb_file, protein_ids):
    """ Fidnthe disordered regions of the proteins in our family to do the analysis on"""
    disordered_regions = {}
    with open(mobidb_file, "r") as file:
        reader = csv.reader(file)
        for row in reader:
            protein_id = row[0]
            # If protein found in family, save its disordered regions
            if protein_id in protein_ids: #
                parsed_list = ast.literal_eval(row[1]) # use literal_eval bc the areas are saved as nested lists [[...]] , and we want to parse it as these lists 
                disordered_regions[protein_id] = [tuple(pair) for pair in parsed_list] # turn each list (i.e. disordered region area) to tuple 
    return disordered_regions

def load_sequences(fasta_file, protein_ids):
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        protein_id = record.id.split("|")[1]
        if protein_id in protein_ids:
            sequences[protein_id] = str(record.seq)
    return sequences

def load_motifs(elm_file, prosite_file):
    motifs = {}
    # Extract needed information from ELM 
    with open(elm_file, 'r') as file:
        for _ in range(5):  # Skip first 5 header lines (see elm file, first 5 are just headers)
            next(file)
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            motif_name = row['ELMIdentifier']
            pattern = row['Regex'].strip()
            motifs[motif_name] = pattern

    # Extract needed information from Prosite
    with open(prosite_file, 'r') as file:
        for line in file:
            # Split on first tab to separate AC from pattern
            parts = line.strip().split('\t') # split on tab character, as entries on one row are like : PS00001 "TAB" N[ACDEFGHIKLMNQRSTVWY][ST][ACDEFGHIKLMNQRSTVWY]
            if len(parts) == 2:
                ac, pattern = parts # split into the two parts (i.e. the name and the pattern)
                motifs[ac] = pattern
    return motifs

def match_motifs(sequences, disordered_regions, motifs):
    # in the sequences found in the family that have a disordered region (or we atleast found one), now look into the disordered regions and see if we find any motifs
    results = {}
    for protein_id, sequence in sequences.items():
        if protein_id in disordered_regions:
            results[protein_id] = []
            for start, end in disordered_regions[protein_id]:
                region_seq = sequence[start - 1:end] # start-1 due to python indexing starting from 0 
                # Now that we are looking at a disordered region in one of the proteins of our family, lets
                # see if we find any motifs
                for motif_name, pattern in motifs.items():
                    # Example values:
                    #pattern = "N[ACDEFGHIKLMNQRSTVWY][ST][ACDEFGHIKLMNQRSTVWY]"
                    #region_seq = "NKSTMNLSTPNQSTV"

                    # The re.finditer() will find ALL occurrences where:
                    # - N matches literally
                    # - followed by ANY ONE of the amino acids in [ACDEFGHIKLMNQRSTVWY]
                    # - followed by either S or T
                    # - followed by ANY ONE of the amino acids in [ACDEFGHIKLMNQRSTVWY]
                    # --> that is what matches does with the regex search
                    matches = re.finditer(r"{}".format(pattern), region_seq)
                    for match in matches:
                        results[protein_id].append({
                            "motif": motif_name,
                            "match": match.group(),
                            "start": start + match.start(),
                            "end": start + match.end()
                        })
    return results

def save_results(results, output_file):
    with open(output_file, "w") as file:
        for protein_id, matches in results.items():
            file.write(f"> {protein_id}\n")
            for match in matches:
                file.write(
                    f"{match['motif']}\t{match['match']}\t{match['start']}-{match['end']}\n"
                )

def main():
    mobidb_file = "mobidb_lite_swissprot.csv"
    fasta_file = "uniprot_sprot.fasta"
    elm_file = "elm_classes.tsv"
    prosite_file = "prosite_preprocessed.txt"
    psiblast_file = "psiblast_parsed.csv"
    hmm_file = "hmmsearch_output.csv"
    output_file = "conserved_motifs_in_disorder.txt"

    protein_ids = load_protein_ids(psiblast_file, hmm_file)
    sequences = load_sequences(fasta_file, protein_ids)
    disordered_regions = load_disordered_regions(mobidb_file, protein_ids)
    motifs = load_motifs(elm_file, prosite_file)

    results = match_motifs(sequences, disordered_regions, motifs)
    save_results(results, output_file)
    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    main()