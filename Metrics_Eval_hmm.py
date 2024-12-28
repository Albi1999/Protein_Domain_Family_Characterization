import pandas as pd
import numpy as np
from math import sqrt

def load_pfam_data(pfam_file):
    """Load and process Pfam domain position data"""
    df = pd.read_csv(pfam_file)
    # Create a set of proteins with Pfam domains
    pfam_proteins = set(df['uniprot_id'])
    # Create dictionary of domain positions for each protein
    pfam_positions = {}
    for _, row in df.iterrows():
        protein = row['uniprot_id']
        start = row['domain_start']
        end = row['domain_end']
        if protein not in pfam_positions:
            pfam_positions[protein] = set()
        # Add all positions in the domain range
        pfam_positions[protein].update(range(start, end + 1))
    return pfam_proteins, pfam_positions

def parse_hmmsearch_output(hmmsearch_file):
    """Parse hmmsearch output to get protein matches and domain positions"""
    hmm_proteins = set()
    hmm_positions = {}
    
    with open(hmmsearch_file, 'r') as f:
        lines = f.readlines()
        
    in_domain_annotation = False
    current_protein = None
    
    for line in lines:
        if line.startswith(">>"):
            in_domain_annotation = False
            # Extract protein ID from header line
            current_protein = line.split()[1].split("|")[1]
            hmm_proteins.add(current_protein)
            if current_protein not in hmm_positions:
                hmm_positions[current_protein] = set()
            
        elif line.strip().startswith("==") and "domain" in line:
            in_domain_annotation = True
            continue
            
        elif in_domain_annotation and line.strip() and not line.startswith("//"):
            # Parse alignment lines to get domain positions
            parts = line.split()
            if len(parts) >= 8 and parts[0] not in ["==", "#"]:
                try:
                    ali_from = int(parts[9])
                    ali_to = int(parts[10])
                    hmm_positions[current_protein].update(range(ali_from, ali_to + 1))
                except (ValueError, IndexError):
                    continue
    
    return hmm_proteins, hmm_positions

def calculate_metrics(tp, fp, tn, fn):
    """Calculate various performance metrics"""
    try:
        precision = tp / (tp + fp) if (tp + fp) > 0 else 0
        recall = tp / (tp + fn) if (tp + fn) > 0 else 0
        f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
        balanced_accuracy = ((tp / (tp + fn)) + (tn / (tn + fp))) / 2 if (tp + fn) > 0 and (tn + fp) > 0 else 0
        
        # Matthews Correlation Coefficient
        numerator = (tp * tn) - (fp * fn)
        denominator = sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)) if (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn) > 0 else 1
        mcc = numerator / denominator
        
        return {
            'Precision': precision,
            'Recall': recall,
            'F1-score': f1,
            'Balanced Accuracy': balanced_accuracy,
            'MCC': mcc
        }
    except ZeroDivisionError:
        return {
            'Precision': 0,
            'Recall': 0,
            'F1-score': 0,
            'Balanced Accuracy': 0,
            'MCC': 0
        }

def evaluate_protein_level(pfam_proteins, hmm_proteins, all_proteins):
    """Evaluate performance at the protein level"""
    tp = len(pfam_proteins & hmm_proteins)
    fp = len(hmm_proteins - pfam_proteins)
    fn = len(pfam_proteins - hmm_proteins)
    tn = len(all_proteins - pfam_proteins - hmm_proteins)
    
    return calculate_metrics(tp, fp, tn, fn)

def evaluate_residue_level(pfam_positions, hmm_positions):
    """Evaluate performance at the residue level"""
    tp = fp = fn = tn = 0
    
    # Get all proteins that appear in either dataset
    all_proteins = set(pfam_positions.keys()) | set(hmm_positions.keys())
    
    for protein in all_proteins:
        pfam_pos = pfam_positions.get(protein, set())
        hmm_pos = hmm_positions.get(protein, set())
        
        # Calculate metrics for this protein
        tp += len(pfam_pos & hmm_pos)
        fp += len(hmm_pos - pfam_pos)
        fn += len(pfam_pos - hmm_pos)
        
        # For TN, we'll consider positions up to the maximum position in either set
        max_pos = max(max(pfam_pos) if pfam_pos else 0, max(hmm_pos) if hmm_pos else 0)
        all_positions = set(range(1, max_pos + 1))
        tn += len(all_positions - pfam_pos - hmm_pos)
    
    return calculate_metrics(tp, fp, tn, fn)

def main():
    # Load Pfam data
    pfam_proteins, pfam_positions = load_pfam_data('pfam_domain_positions.csv')
    
    # Parse HMMsearch output
    hmm_proteins, hmm_positions = parse_hmmsearch_output('hmmsearch_output.txt')
    
    # Create set of all proteins (universe)
    all_proteins = pfam_proteins | hmm_proteins

    print("HMMSEARCH vs Pfam Evaluation Results")
    print("-" * 40)
    
    # Evaluate at protein level
    print("Protein-level evaluation:")
    protein_metrics = evaluate_protein_level(pfam_proteins, hmm_proteins, all_proteins)
    for metric, value in protein_metrics.items():
        print(f"{metric}: {value:.3f}")
    
    print("\nResidue-level evaluation:")
    residue_metrics = evaluate_residue_level(pfam_positions, hmm_positions)
    for metric, value in residue_metrics.items():
        print(f"{metric}: {value:.3f}")

if __name__ == "__main__":
    main()