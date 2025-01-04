import re
import csv
from Bio import SeqIO
import ast

def load_disordered_regions(mobidb_file, pfam_ids):
    disordered_regions = {}
    with open(mobidb_file, "r") as file:
        reader = csv.reader(file)
        for row in reader:
            protein_id = row[0]
            if protein_id in pfam_ids:
                parsed_list = ast.literal_eval(row[1])
                disordered_regions[protein_id] = [tuple(pair) for pair in parsed_list]
    return disordered_regions

def load_sequences(fasta_file):
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        protein_id = record.id.split("|")[1]
        sequences[protein_id] = str(record.seq)
    return sequences

def match_motifs(sequences, disordered_regions, motifs):
    results = {}
    for protein_id, sequence in sequences.items():
        if protein_id in disordered_regions:
            results[protein_id] = []
            for start, end in disordered_regions[protein_id]:
                region_seq = sequence[start - 1:end]
                for motif_name, pattern in motifs.items():

                    matches = re.finditer(r"{}".format(pattern), region_seq)
                    for match in matches:
                        results[protein_id].append({
                            "motif": motif_name,
                            "match": match.group(),
                            "start": start + match.start(),
                            "end": start + match.end()
                        })
    return results

def load_motifs(elm_file, prosite_file):
    motifs = {}
    with open(elm_file, 'r') as file:
        for _ in range(5):  # Skip first 5 header lines
            next(file)
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            motif_name = row['ELMIdentifier']
            pattern = row['Regex'].strip()
            motifs[motif_name] = pattern

    with open(prosite_file, 'r') as file:
        for line in file:
            # Split on first tab to separate AC from pattern
            parts = line.strip().split('\t')
            if len(parts) == 2:
                ac, pattern = parts
                motifs[ac] = pattern
    return motifs


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
    fasta_file = "uniprotkb_xref_pfam_pf00151_AND_reviewe_2024_12_31.fasta"
    elm_file = "elm_classes.tsv"
    prosite_file = "prosite_preprocessed.txt"
    output_file = "conserved_motifs_in_disorder.txt"

    sequences = load_sequences(fasta_file)
    pfam_ids = set(sequences.keys())
    disordered_regions = load_disordered_regions(mobidb_file, pfam_ids)
    motifs = load_motifs(elm_file, prosite_file)

    results = match_motifs(sequences, disordered_regions, motifs)
    save_results(results, output_file)
    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    main()