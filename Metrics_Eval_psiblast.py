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

def parse_psiblast_output(psiblast_file):
    """Parse PSI-BLAST output to get protein matches and alignment positions"""
    # Column names for the PSI-BLAST output
    columns = ['qseqid', 'sseqid', 'qstart', 'qend', 'sstart', 'send', 'pident', 'evalue']
    
    # Read the PSI-BLAST output
    blast_df = pd.read_csv(psiblast_file, sep='\t', names=columns)
    
    # Extract protein IDs from sseqid column (format: sp|XXXXX|NAME)
    blast_df['protein_id'] = blast_df['sseqid'].apply(lambda x: x.split('|')[1])
    
    # Create set of proteins found by PSI-BLAST
    blast_proteins = set(blast_df['protein_id'])
    
    # Create dictionary of matched positions for each protein
    blast_positions = {}
    for _, row in blast_df.iterrows():
        protein = row['protein_id']
        start = int(row['sstart'])
        end = int(row['send'])
        
        if protein not in blast_positions:
            blast_positions[protein] = set()
        
        # Add all positions in the alignment range
        blast_positions[protein].update(range(start, end + 1))
    
    return blast_proteins, blast_positions

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

def evaluate_protein_level(pfam_proteins, blast_proteins, all_proteins):
    """Evaluate performance at the protein level"""
    tp = len(pfam_proteins & blast_proteins)
    fp = len(blast_proteins - pfam_proteins)
    fn = len(pfam_proteins - blast_proteins)
    tn = len(all_proteins - pfam_proteins - blast_proteins)
    
    return calculate_metrics(tp, fp, tn, fn)

def evaluate_residue_level(pfam_positions, blast_positions):
    """Evaluate performance at the residue level"""
    tp = fp = fn = tn = 0
    
    # Get all proteins that appear in either dataset
    all_proteins = set(pfam_positions.keys()) | set(blast_positions.keys())
    
    for protein in all_proteins:
        pfam_pos = pfam_positions.get(protein, set())
        blast_pos = blast_positions.get(protein, set())
        
        # Calculate metrics for this protein
        tp += len(pfam_pos & blast_pos)
        fp += len(blast_pos - pfam_pos)
        fn += len(pfam_pos - blast_pos)
        
        # For TN, we'll consider positions up to the maximum position in either set
        max_pos = max(max(pfam_pos) if pfam_pos else 0, max(blast_pos) if blast_pos else 0)
        all_positions = set(range(1, max_pos + 1))
        tn += len(all_positions - pfam_pos - blast_pos)
    
    return calculate_metrics(tp, fp, tn, fn)

def main():
    # Load Pfam data
    pfam_proteins, pfam_positions = load_pfam_data('pfam_domain_positions.csv')
    
    # Parse PSI-BLAST output
    blast_proteins, blast_positions = parse_psiblast_output('psiblast_search_output.txt')
    
    # Create set of all proteins (universe)
    all_proteins = pfam_proteins | blast_proteins
    
    print("PSI-BLAST vs Pfam Evaluation Results")
    print("-" * 40)
    
    # Evaluate at protein level
    print("\nProtein-level evaluation:")
    protein_metrics = evaluate_protein_level(pfam_proteins, blast_proteins, all_proteins)
    for metric, value in protein_metrics.items():
        print(f"{metric}: {value:.3f}")
    
    print("\nResidue-level evaluation:")
    residue_metrics = evaluate_residue_level(pfam_positions, blast_positions)
    for metric, value in residue_metrics.items():
        print(f"{metric}: {value:.3f}")

if __name__ == "__main__":
    main()