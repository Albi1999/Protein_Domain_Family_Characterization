from Bio import AlignIO
from collections import Counter
import pandas as pd
from scipy.stats import entropy
import math
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
import sys

class ConservationAnalyzer:
    def __init__(self, alignment_file):
        """
        Initialize with an alignment file
            alignment_file (str): Path to the alignment file
        """
        self.alignment = AlignIO.read(alignment_file, 'fasta')
        self.num_sequences = len(self.alignment)
        self.alignment_length = self.alignment.get_alignment_length()
        
    def get_column(self, pos):
        """Extract a column from the alignment"""
        return [record.seq[pos] for record in self.alignment]
    
    def calculate_gap_frequency(self, pos):
        """Calculate frequency of gaps in a column"""
        column = self.get_column(pos)
        return column.count('-') / len(column)
    
    def calculate_amino_acid_frequencies(self, pos):
        """Calculate frequencies of each amino acid in a column"""
        column = self.get_column(pos)
        total = len(column) - column.count('-')  # Don't count gaps, such that when we calculate conservation scores the gaps don't mess it up 
        if total == 0:
            return {}
        
        counts = Counter(aa for aa in column if aa != '-')
        return {aa: count/total for aa, count in counts.items()}
    
    def calculate_conservation_score(self, pos):
        """
        Calculate conservation score based on frequency of most common amino acid
        Ignores gaps in calculation
        """
        freqs = self.calculate_amino_acid_frequencies(pos)
        if not freqs:
            return 0
        return max(freqs.values())
    
    def calculate_entropy(self, pos):
        """
        Calculate Shannon entropy for a column
        Lower entropy means higher conservation
        """
        freqs = self.calculate_amino_acid_frequencies(pos)
        if not freqs:
            return float('inf')  
        
        return -sum(p * math.log2(p) for p in freqs.values())
    
    def get_amino_acid_groups(self):
        """Define groups of similar amino acids 
           Based on : https://en.wikipedia.org/wiki/Conservative_replacement#:~:text=There%20are%2020%20naturally%20occurring,both%20small%2C%20negatively%20charged%20residues.
        """
        return {
            'aliphatic': set('GAVLI'),
            'hydroxyl': set('SCUTM'),
            'cyclic': set('P'),
            'aromatic': set('FYW'),
            'basic': set('HKR'),
            'acidic': set('DENQ')
        }
    
    def calculate_group_conservation(self, pos):
        """
        Calculate conservation considering amino acid groups
        Basically the same as calculate_conversation_score, just that it calculates based on the groups, not single amino acids !
        """
        column = self.get_column(pos)
        groups = self.get_amino_acid_groups()
        
        # Assign each amino acid to its group
        aa_to_group = {}
        for group_name, aas in groups.items():
            for aa in aas:
                aa_to_group[aa] = group_name
        
        # Count group occurrences
        group_counts = Counter(aa_to_group.get(aa, 'other') 
                             for aa in column if aa != '-')
        
        if not group_counts:
            return 0
            
        return max(group_counts.values()) / sum(group_counts.values())





    # TODO : I took very strict values now such that the number of residues per sequence is below 100 (right now we have length 77) ; the PSSM creation with 
    # much higher length did not work, but maybe we should write an email and ask ; nevertheless, we can first try some evaluation based on that PSSM and see our scores

    # TODO : diff gap_thresh/conservation threshold for different number of columns in output (OPTIMIZE)
    def analyze_columns(self, gap_threshold=0.37, conservation_threshold=0.9):
        """
        Analyze all columns and return comprehensive metrics
        Returns DataFrame with various conservation metrics for each position
        """
        data = []
        
        for i in range(self.alignment_length):
            gap_freq = self.calculate_gap_frequency(i)
            cons_score = self.calculate_conservation_score(i)
            info_content = self.calculate_entropy(i)
            group_cons = self.calculate_group_conservation(i)
            
            data.append({
                'position': i + 1,
                'gap_frequency': gap_freq,
                'single_conservation': cons_score,
                'entropy': info_content,
                'group_conservation': group_cons,
                # Here we should look possibly for better ideas
                # Check gap frequency not too high (i.e. not nearly all elements in the columns gaps (-))
                # Check that the group conservation is high enough (i.e. the amino acids are not too different
                # ; right now we do with groups and not single amino acid sequence since I'd say the groups
                # are more representative (if we do single amino acids, we'd delete more stuff))
                'suggested_remove': (gap_freq > gap_threshold or       
                                   group_cons < conservation_threshold) # TODO : OPTIMIZE WHEN TO REMOVE
            })
        
        return pd.DataFrame(data)


def remove_columns_from_alignment(input_file, output_file, columns_to_remove, format="fasta"):
    """
    Remove specified columns from a multiple sequence alignment and save to new file
    
    Args:
        input_file (str): Path to input alignment file
        output_file (str): Path where to save trimmed alignment
        columns_to_remove (list): List of column indices to remove (0-based)
        format (str): File format (default: "fasta")
    """
    # Read the alignment
    alignment = AlignIO.read(input_file, format)
    
    # Sort columns to remove in descending order
    # (so removing them doesn't affect the indices of remaining columns)
    columns_to_remove = sorted(columns_to_remove, reverse=True)
    
    # Create new alignment records
    new_records = []
    
    # Process each sequence
    for record in alignment:
        # Convert sequence to list for easier manipulation
        seq_list = list(record.seq)
        
        # Remove specified columns
        for col in columns_to_remove:
            del seq_list[col]
        
        # Create new sequence record
        new_seq = Seq(''.join(seq_list)) # Join the list element to a string again (i.e. after removal of amino acids out of sequence represented as list, turn into one string again) and turn into Seq object
        new_record = SeqRecord(new_seq,
                            id=record.id,
                            name=record.name,
                            description=record.description)
        new_records.append(new_record)
    
    # Create new alignment
    # TODO : Maybe we have to add some variables here (i.e. how to do the MSA)!
    new_alignment = MultipleSeqAlignment(new_records)
    
    # Write to file
    AlignIO.write(new_alignment, output_file, format)
    
    return new_alignment




    


# Example usage:
if __name__ == "__main__":
    # Initialize analyzer 
    analyzer = ConservationAnalyzer("clustal_rows_removed_100threshold.fa")
    
    # Get comprehensive analysis
    analysis = analyzer.analyze_columns()
   # analysis_2 = analyzer.analyze_rows()
    
    # Print summary statistics
    print("\nAlignment Summary:")
    print(f"Number of sequences: {analyzer.num_sequences}")
    print(f"Alignment length: {analyzer.alignment_length}")


    # Print number of True/False
    counts = analysis['suggested_remove'].value_counts()

    counts_true = counts[True]  # To be removed
    counts_false = counts[False] # To be kept

    print(f"With the current removal tactic, we would remove {(counts_true / (counts_true + counts_false)):.2f} percent of columns ; we keep {counts_false} of {counts_false + counts_true} columns")
    

    # Save detailed analysis to CSV
    analysis.to_csv("conservation_analysis.csv", index=False)


    # Get indices of columns marked for removal
    columns_to_remove = analysis[analysis['suggested_remove']]['position'].values.tolist()
    # Convert to 0-based indices (if positions were 1-based)
    columns_to_remove = [x-1 for x in columns_to_remove]
    
    # Remove columns and save new alignment
    new_alignment = remove_columns_from_alignment(
        "clustal_rows_removed_100threshold.fa",
        "trimmed_alignment.fasta",
        columns_to_remove
    )


        


    print(f"Original alignment length: {analyzer.alignment_length}")
    print(f"Number of columns removed: {len(columns_to_remove)}")
    print(f"New alignment length: {new_alignment.get_alignment_length()}")


