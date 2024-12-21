from Bio import AlignIO
from collections import Counter
import pandas as pd
from scipy.stats import entropy
import math


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
    
    def analyze_columns(self, gap_threshold=0.5, conservation_threshold=0.5):
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
                'conservation_score': cons_score,
                'entropy': info_content,
                'group_conservation': group_cons,
                # Here we should look possibly for better ideas
                'suggested_remove': (gap_freq > gap_threshold or       
                                   group_cons < conservation_threshold)
            })
        
        return pd.DataFrame(data)
    
    """
    def suggest_columns_to_remove(self, gap_threshold=0.5, 
                                conservation_threshold=0.3,
                                information_threshold=0.5):
      
        to_remove = []
        
        for i in range(self.alignment_length):
            gap_freq = self.calculate_gap_frequency(i)
            cons_score = self.calculate_conservation_score(i)
            info_content = self.calculate_entropy(i)
            
            # Suggest removal if:
            # 1. Too many gaps OR
            # 2. Low conservation AND low information content
            if (gap_freq > gap_threshold or 
                (cons_score < conservation_threshold and 
                 info_content < information_threshold)):
                to_remove.append(i)
                
        return to_remove
    """

# Example usage:
if __name__ == "__main__":
    # Initialize analyzer 
    analyzer = ConservationAnalyzer("clustalo-I20241221-155509-0013-19108079-p1m.fasta")
    
    # Get comprehensive analysis
    analysis = analyzer.analyze_columns()
    
    # Print summary statistics
    print("\nAlignment Summary:")
    print(f"Number of sequences: {analyzer.num_sequences}")
    print(f"Alignment length: {analyzer.alignment_length}")
    
    # Get suggested columns to remove
  #  to_remove = analyzer.suggest_columns_to_remove()
  #  print(f"\nSuggested columns to remove: {len(to_remove)}")

    # Print number of True/False
    counts = analysis['suggested_remove'].value_counts()

    counts_true = counts[True]  # To be removed
    counts_false = counts[False] # To be kept

    print(f"With the current removal tactic, we would remove {(counts_true / (counts_true + counts_false)):.2f} percent of columns")
    

    
    # Save detailed analysis to CSV
    analysis.to_csv("conservation_analysis.csv", index=False)

    # TODO : ADD REMOVAL LOGIC OF ROWS