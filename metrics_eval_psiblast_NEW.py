"""Based on the .csv files that closely resemble each other, do the calculations"""

import pandas as pd
import numpy as np
from sklearn.metrics import precision_score, recall_score, f1_score, balanced_accuracy_score, matthews_corrcoef



def create_residue_vectors(psiblast_df, pfam_df, protein_id):
    """
    For a single protein (based on its protein_id), create residue vectors
    """
    # Get domain positions from both predictions
    psi_matches = psiblast_df[psiblast_df['uniprot_id'] == protein_id]
    pfam_matches = pfam_df[pfam_df['uniprot_id'] == protein_id]
    
    # Check if for the protein at hand, we even find it in both of the .csv's 
    if len(psi_matches) == 0 or len(pfam_matches) == 0:
        return None
    
    # Get the max of the max lengths for each hit for the current protein
    # Notice that we can have multiple hits per protein (i.e. multiple alignments that were found)
    # So we have to account for all of them
    max_length = max(
        psi_matches['domain_end'].max(),
        pfam_matches['domain_end'].max()
    )
    
    # With that, we can create vectors of the same size
    # Create binary vectors for each position
    true_positions = np.zeros(max_length)
    pred_positions = np.zeros(max_length)
    

    # Now, iterate through all the found alignments (for the current protein)

    # Fill in Pfam (true) positions
    for _, row in pfam_matches.iterrows():
        start = row['domain_start'] - 1  # Convert to 0-based indexing
        end = row['domain_end']
        true_positions[start:end] = 1
        
    # Fill in PSIBLAST (predicted) positions
    for _, row in psi_matches.iterrows():
        start = row['domain_start'] - 1  # Convert to 0-based indexing
        end = row['domain_end']
        pred_positions[start:end] = 1
    
    return true_positions, pred_positions

def evaluate_model_protein_level(psiblast_file, pfam_file, only_found=False):
    """
    Evaluate PSIBLAST model performance against Pfam annotations
    
    Parameters:
    - psiblast_file: Path to PSIBLAST results CSV
    - pfam_file: Path to Pfam ground truth CSV
    - only_found: If True, only evaluate proteins found by PSIBLAST
    """
    # Step 1: Load both CSV files
    psiblast_df = pd.read_csv(psiblast_file)
    pfam_df = pd.read_csv(pfam_file)
    
    # Step 2: Get unique list of proteins from both files
    psiblast_proteins = set(psiblast_df['uniprot_id'])
    pfam_proteins = set(pfam_df['uniprot_id'])
    # Use union operator to get all the proteins in total
    all_proteins = psiblast_proteins.union(pfam_proteins)
    
    if only_found:
        # Only consider proteins that PSIBLAST found
        all_proteins = psiblast_proteins
        print("\nEvaluating only PSIBLAST-found proteins:")
    else:
        # Consider all proteins from both sets
        all_proteins = psiblast_proteins.union(pfam_proteins)
        print("\nEvaluating all proteins:")
    
    print(f"Number of proteins predicted by PSIBLAST: {len(psiblast_proteins)}")
    print(f"Number of proteins in Pfam ground truth: {len(pfam_proteins)}")
    print(f"Number of proteins being evaluated: {len(all_proteins)}")
    

    print("\n=== Protein-Level Evaluation ===")
    # Step 3: Create binary vectors for true and predicted labels
    y_true = []  # Ground truth from Pfam
    y_pred = []  # Predictions from PSIBLAST
    
    for protein in all_proteins:
        y_true.append(1 if protein in pfam_proteins else 0)
        y_pred.append(1 if protein in psiblast_proteins else 0)

    # So we have something like
    # y_true 0 0 1 0 1 ...
    # y_pred 0 1 1 0 1 ...
    
    # Step 4: Calculate performance metrics
    protein_results = {
        'Precision': precision_score(y_true, y_pred),
        'Recall': recall_score(y_true, y_pred),
        'F-score': f1_score(y_true, y_pred),
        'Balanced Accuracy': balanced_accuracy_score(y_true, y_pred),
        'MCC': matthews_corrcoef(y_true, y_pred)
    }


    # Step 5: Calculate confusion matrix components
    tp = sum(1 for t, p in zip(y_true, y_pred) if t == 1 and p == 1)
    fp = sum(1 for t, p in zip(y_true, y_pred) if t == 0 and p == 1)
    fn = sum(1 for t, p in zip(y_true, y_pred) if t == 1 and p == 0)
    tn = sum(1 for t, p in zip(y_true, y_pred) if t == 0 and p == 0)
    
    # Print detailed results
    print("\nProtein-Level Confusion Matrix:")
    print(f"True Positives: {tp}")
    print(f"False Positives: {fp}")
    print(f"False Negatives: {fn}")
    print(f"True Negatives: {tn}")


    print("\nProtein-Level Metrics:")
    for metric, value in protein_results.items():
        print(f"{metric}: {value:.4f}")
    
    # Residue-level evaluation
    print("\n=== Residue-Level Evaluation ===")
    # Only evaluate residues for proteins found in both sets 
    common_proteins = psiblast_proteins.intersection(pfam_proteins)
    print(f"Number of proteins for residue-level evaluation: {len(common_proteins)}")
    
    # Collect all residue-level predictions
    all_true_residues = []
    all_pred_residues = []
    
    for protein in common_proteins:
        result = create_residue_vectors(psiblast_df, pfam_df, protein)
        if result is not None:
            true_pos, pred_pos = result
            all_true_residues.extend(true_pos)
            all_pred_residues.extend(pred_pos)
    
    # Calculate residue-level metrics
    residue_results = {
        'Precision': precision_score(all_true_residues, all_pred_residues),
        'Recall': recall_score(all_true_residues, all_pred_residues),
        'F-score': f1_score(all_true_residues, all_pred_residues),
        'Balanced Accuracy': balanced_accuracy_score(all_true_residues, all_pred_residues),
        'MCC': matthews_corrcoef(all_true_residues, all_pred_residues)
    }
    
    # Calculate residue-level confusion matrix
    tp = sum(1 for t, p in zip(all_true_residues, all_pred_residues) if t == 1 and p == 1)
    fp = sum(1 for t, p in zip(all_true_residues, all_pred_residues) if t == 0 and p == 1)
    fn = sum(1 for t, p in zip(all_true_residues, all_pred_residues) if t == 1 and p == 0)
    tn = sum(1 for t, p in zip(all_true_residues, all_pred_residues) if t == 0 and p == 0)
    
    print("\nResidue-Level Confusion Matrix:")
    print(f"True Positives: {tp}")
    print(f"False Positives: {fp}")
    print(f"False Negatives: {fn}")
    print(f"True Negatives: {tn}")
    
    print("\nResidue-Level Metrics:")
    for metric, value in residue_results.items():
        print(f"{metric}: {value:.4f}")
    
    return protein_results, residue_results


psiblast_file = 'psiblast_parsed.csv'
pfam_file = 'pfam_domain_positions.csv'

#print("===============================")
#print("Evaluation on all proteins:")
#results_all = evaluate_model_protein_level(psiblast_file, pfam_file, only_found=False)

print("\n===============================")
print("Evaluation only on PSIBLAST-found proteins:")
results_found = evaluate_model_protein_level(psiblast_file, pfam_file, only_found=True)