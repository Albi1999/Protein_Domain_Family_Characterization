"""Based on the .csv files that closely resemble each other, do the calculations"""

import pandas as pd
import numpy as np
from sklearn.metrics import precision_score, recall_score, f1_score, balanced_accuracy_score, matthews_corrcoef



def create_residue_vectors(pred_df, pfam_df, protein_id):
    """
    For a single protein (based on its protein_id), create residue vectors
    """

    # Note that since the ground truth (i.e. the PFAM domains) are all just a single sequence hit, we only include the strongest
    # (lowest e-value) hit from the HMM search (also we noticed that basically all 2nd, 3rd domain hits have much higher e-values)
 
    # Get domain positions from both predictions
    pred_matches = pred_df[pred_df['uniprot_id'] == protein_id]
    pfam_matches = pfam_df[pfam_df['uniprot_id'] == protein_id]
    
    # Check if for the protein at hand, we even find it in both of the .csv's 
    if len(pred_matches) == 0 or len(pfam_matches) == 0:
        return None
    

    assert (len(pred_df == 1) , "not length 1 ")
    # Get the max of the max lengths for each hit for the current protein
    # Notice that we can have multiple hits per protein (i.e. multiple alignments that were found)
    # So we have to account for all of them
    max_length = max(
        pred_matches['domain_end'].max(),
        pfam_matches['domain_end'].max()
    )
    
    # With that, we can create vectors of the same size
    # Create binary vectors for each position
    true_positions = np.zeros(int(max_length))
    pred_positions = np.zeros(int(max_length))
    

    # Now, iterate through all the found alignments (for the current protein)

    # Fill in Pfam (true) positions
    for _, row in pfam_matches.iterrows():
        start = row['domain_start'] - 1  # Convert to 0-based indexing
        end = row['domain_end']
        true_positions[start:end] = 1
        
    # Fill in PSSM/HMM (predicted) positions
    for _, row in pred_matches.iterrows():
        start = int(row['domain_start'] - 1)  # Convert to 0-based indexing
        end = int(row['domain_end'])
        pred_positions[start:end] = 1
    
    return true_positions, pred_positions




def evaluate_model(psiblast_file, hmm_file, pfam_file, only_found=False, e_threshold = 0.0001):
    """
    Evaluate PSIBLAST model performance against Pfam annotations
    
    Parameters:
    - psiblast_file: Path to PSIBLAST results CSV
    - hmm_file : Path to HMM results CSV
    - pfam_file: Path to Pfam ground truth CSV
    - only_found: If True, only evaluate proteins found by PSIBLAST
    """
    # Step 1: Load both CSV files
    psiblast_df = pd.read_csv(psiblast_file)
    hmm_df = pd.read_csv(hmm_file)
    pfam_df = pd.read_csv(pfam_file)

    # HMM finds a lot of hits, a lot with extremely high e-values : Take only the ones that are above some threshold (score for now, later e-value)
    # We filter based on the first domain hit (TODO : i.e. the best one ? check if first found domain always strongest)
    filtered_hmm_proteins = hmm_df[hmm_df['E-value'] <= e_threshold]['uniprot_id']
    
    # Step 2: Get unique list of proteins from both files
    psiblast_proteins = set(psiblast_df['uniprot_id'])
    hmm_proteins = set(filtered_hmm_proteins)
    pfam_proteins = set(pfam_df['uniprot_id'])

    
    if only_found:
        # Only consider proteins that PSIBLAST/HMM found
        all_proteins_psiblast = psiblast_proteins
        all_proteins_hmm = hmm_proteins
        print("\nEvaluating only PSIBLAST-found proteins:")
    else:
        # Consider all proteins from both sets
        all_proteins_psiblast = psiblast_proteins.union(pfam_proteins)
        all_proteins_hmm = hmm_proteins.union(pfam_proteins)
        print("\nEvaluating all proteins:")
    
    print(f"Number of proteins predicted by PSIBLAST: {len(psiblast_proteins)}")
    print(f"Number of proteins predicted by HMM: {len(hmm_proteins)}")
    print(f"Number of proteins in Pfam ground truth: {len(pfam_proteins)}")
    print(f"Number of proteins being evaluated for PSIBLAST: {len(all_proteins_psiblast)}")
    print(f"Number of proteins being evaluated for HMM: {len(all_proteins_hmm)}")
    

    print("\n=== Protein-Level Evaluation ===")
    # Step 3: Create binary vectors for true and predicted labels
    y_true_psiblast = []  # Ground truth from Pfam
    y_pred_psiblast = []  # Predictions from PSIBLAST
    y_true_hmm = []
    y_pred_hmm = []


    
    for protein in all_proteins_psiblast:
        y_true_psiblast.append(1 if protein in pfam_proteins else 0)
        y_pred_psiblast.append(1 if protein in psiblast_proteins else 0)


    for protein in all_proteins_hmm:
        y_true_hmm.append(1 if protein in pfam_proteins else 0)
        y_pred_hmm.append(1 if protein in hmm_proteins else 0)

    # So we have something like
    # y_true_psiblast 0 0 1 0 1 ...
    # y_pred_psiblast 0 1 1 0 1 ...
    
    # Step 4: Calculate performance metrics
    protein_results_psiblast = {
        'Precision': precision_score(y_true_psiblast, y_pred_psiblast),
        'Recall': recall_score(y_true_psiblast, y_pred_psiblast),
        'F-score': f1_score(y_true_psiblast, y_pred_psiblast),
        'Balanced Accuracy': balanced_accuracy_score(y_true_psiblast, y_pred_psiblast),
        'MCC': matthews_corrcoef(y_true_psiblast, y_pred_psiblast)
    }


    protein_results_hmm = {
        'Precision': precision_score(y_true_hmm, y_pred_hmm),
        'Recall': recall_score(y_true_hmm, y_pred_hmm),
        'F-score': f1_score(y_true_hmm, y_pred_hmm),
        'Balanced Accuracy': balanced_accuracy_score(y_true_hmm, y_pred_hmm),
        'MCC': matthews_corrcoef(y_true_hmm, y_pred_hmm)
    }


    # Step 5: Calculate confusion matrix components
    tp_psiblast = sum(1 for t, p in zip(y_true_psiblast, y_pred_psiblast) if t == 1 and p == 1)
    fp_psiblast = sum(1 for t, p in zip(y_true_psiblast, y_pred_psiblast) if t == 0 and p == 1)
    fn_psiblast = sum(1 for t, p in zip(y_true_psiblast, y_pred_psiblast) if t == 1 and p == 0)
    tn_psiblast = sum(1 for t, p in zip(y_true_psiblast, y_pred_psiblast) if t == 0 and p == 0)


    tp_hmm = sum(1 for t, p in zip(y_true_hmm, y_pred_hmm) if t == 1 and p == 1)
    fp_hmm = sum(1 for t, p in zip(y_true_hmm, y_pred_hmm) if t == 0 and p == 1)
    fn_hmm = sum(1 for t, p in zip(y_true_hmm, y_pred_hmm) if t == 1 and p == 0)
    tn_hmm = sum(1 for t, p in zip(y_true_hmm, y_pred_hmm) if t == 0 and p == 0)
    
    # Print detailed results
    print("\nProtein-Level Confusion Matrix PSIBLAST:")
    print(f"True Positives: {tp_psiblast}")
    print(f"False Positives: {fp_psiblast}")
    print(f"False Negatives: {fn_psiblast}")
    print(f"True Negatives: {tn_psiblast}")



    print("\nProtein-Level Confusion Matrix HMM:")
    print(f"True Positives: {tp_hmm}")
    print(f"False Positives: {fp_hmm}")
    print(f"False Negatives: {fn_hmm}")
    print(f"True Negatives: {tn_hmm}")


    print("\nProtein-Level Metrics PSIBLAST:")
    for metric, value in protein_results_psiblast.items():
        print(f"{metric}: {value:.4f}")


    print("\nProtein-Level Metrics HMM:")
    for metric, value in protein_results_hmm.items():
        print(f"{metric}: {value:.4f}")



    # Residue-level evaluation
   
    print("\n=== Residue-Level Evaluation ===")
    # Only evaluate residues for proteins found in both sets 
    common_proteins_psiblast = psiblast_proteins.intersection(pfam_proteins)
    common_proteins_hmm = hmm_proteins.intersection(pfam_proteins)
    print(f"Number of proteins for residue-level evaluation PSIBLAST: {len(common_proteins_psiblast)}")
    print(f"Number of proteins for residue-level evaluation HMM: {len(common_proteins_hmm)}")
    
    # Collect all residue-level predictions
    all_true_residues_psiblast = []
    all_pred_residues_psiblast = []

    all_true_residues_hmm = []
    all_pred_residues_hmm = []

    
    for protein in common_proteins_psiblast:
        result = create_residue_vectors(psiblast_df, pfam_df, protein)
        if result is not None:
            true_pos, pred_pos = result
            all_true_residues_psiblast.extend(true_pos)
            all_pred_residues_psiblast.extend(pred_pos)


    for protein in common_proteins_hmm:
        result = create_residue_vectors(hmm_df, pfam_df, protein)
        if result is not None:
            true_pos, pred_pos = result
            all_true_residues_hmm.extend(true_pos)
            all_pred_residues_hmm.extend(pred_pos)
    
    # Calculate residue-level metrics
    residue_results_psiblast = {
        'Precision': precision_score(all_true_residues_psiblast, all_pred_residues_psiblast),
        'Recall': recall_score(all_true_residues_psiblast, all_pred_residues_psiblast),
        'F-score': f1_score(all_true_residues_psiblast, all_pred_residues_psiblast),
        'Balanced Accuracy': balanced_accuracy_score(all_true_residues_psiblast, all_pred_residues_psiblast),
        'MCC': matthews_corrcoef(all_true_residues_psiblast, all_pred_residues_psiblast)
    }
    
    # Calculate residue-level confusion matrix
    tp = sum(1 for t, p in zip(all_true_residues_psiblast, all_pred_residues_psiblast) if t == 1 and p == 1)
    fp = sum(1 for t, p in zip(all_true_residues_psiblast, all_pred_residues_psiblast) if t == 0 and p == 1)
    fn = sum(1 for t, p in zip(all_true_residues_psiblast, all_pred_residues_psiblast) if t == 1 and p == 0)
    tn = sum(1 for t, p in zip(all_true_residues_psiblast, all_pred_residues_psiblast) if t == 0 and p == 0)
    
    print("\nResidue-Level Confusion Matrix PSIBLAST:")
    print(f"True Positives: {tp}")
    print(f"False Positives: {fp}")
    print(f"False Negatives: {fn}")
    print(f"True Negatives: {tn}")
    
    print("\nResidue-Level Metrics PSIBLAST:")
    for metric, value in residue_results_psiblast.items():
        print(f"{metric}: {value:.4f}")



    # Calculate residue-level metrics
    residue_results_hmm = {
        'Precision': precision_score(all_true_residues_hmm, all_pred_residues_hmm),
        'Recall': recall_score(all_true_residues_hmm, all_pred_residues_hmm),
        'F-score': f1_score(all_true_residues_hmm, all_pred_residues_hmm),
        'Balanced Accuracy': balanced_accuracy_score(all_true_residues_hmm, all_pred_residues_hmm),
        'MCC': matthews_corrcoef(all_true_residues_hmm, all_pred_residues_hmm)
    }
    
    # Calculate residue-level confusion matrix
    tp = sum(1 for t, p in zip(all_true_residues_hmm, all_pred_residues_hmm) if t == 1 and p == 1)
    fp = sum(1 for t, p in zip(all_true_residues_hmm, all_pred_residues_hmm) if t == 0 and p == 1)
    fn = sum(1 for t, p in zip(all_true_residues_hmm, all_pred_residues_hmm) if t == 1 and p == 0)
    tn = sum(1 for t, p in zip(all_true_residues_hmm, all_pred_residues_hmm) if t == 0 and p == 0)


    print("\nResidue-Level Confusion Matrix HMM:")
    print(f"True Positives: {tp}")
    print(f"False Positives: {fp}")
    print(f"False Negatives: {fn}")
    print(f"True Negatives: {tn}")
    
    print("\nResidue-Level Metrics HMM:")
    for metric, value in residue_results_hmm.items():
        print(f"{metric}: {value:.4f}")
    



psiblast_file = 'psiblast_parsed.csv'
hmm_file = 'hmmsearch_output.csv'
pfam_file = 'pfam_domain_positions.csv'

#print("===============================")
#print("Evaluation on all proteins:")
#results_all = evaluate_model_protein_level(psiblast_file, pfam_file, only_found=False)

print("\n===============================")
print("Evaluation only on found proteins in both PSSM/HMM:")
evaluate_model(psiblast_file,hmm_file, pfam_file, only_found=False, e_threshold= 0.001)