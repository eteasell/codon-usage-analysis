'''
Static functions for accessing data files.
'''

import pandas as pd
from Bio import SeqIO

def load_elongation_rates():
    path = '../data/daoduc_song_2018_rates.csv'
    return pd.read_csv(path, index_col='Codon')
    

def load_rscu_data(path):
    return pd.read_csv(path, sep='\t', index_col='codon')
    

def save_rscu(path, window_rscu, full_rscu):
    
    df = pd.DataFrame({
        "window": window_rscu,
        "complement": full_rscu
    })
    df.index.name = 'codon'
    df.to_csv(path, index=True, sep='\t')
    print(f"Saved RSCU data to {path}")
    
    
def load_fasta_records(path: str, genes: list[str] | None = None):
    records = list(SeqIO.parse(path, 'fasta'))
    
    if genes is not None:
        filtered_records = []
        for record in records:
            if record.id in genes:
                filtered_records.append(record)
        return filtered_records
    else:
        return records
                
          
########## For reading GFF3 outputs from TMHMM ############          
      
def parse_gff3_tmr_file(filepath):
    """
    Parse a GFF3 file containing transmembrane region predictions.
    
    Args:
        filepath (str): Path to the GFF3 file
        
    Returns:
        list: List of dictionaries, each containing data for one sequence:
            - 'sequence_id': sequence identifier
            - 'length': sequence length
            - 'tmr_count': number of predicted transmembrane regions
            - 'features': list of feature dictionaries with keys:
                - 'seqid': sequence ID
                - 'source': source (always 'inside', 'outside', or 'TMhelix')
                - 'type': feature type (same as source in this format)
                - 'start': start position (1-based)
                - 'end': end position (1-based)
                - 'score': score (empty in this format)
                - 'strand': strand (empty in this format)
                - 'phase': phase (empty in this format)
                - 'attributes': attributes (empty in this format)
    """
    sequences = []
    current_sequence = None
    
    with open(filepath, 'r') as file:
        for line in file:
            line = line.strip()
            
            # Skip empty lines
            if not line:
                continue
                
            # Handle version header
            if line.startswith('##gff-version'):
                continue
                
            # Handle sequence separator
            if line == '//':
                if current_sequence:
                    sequences.append(current_sequence)
                current_sequence = None
                continue
                
            # Handle comment lines with sequence metadata
            if line.startswith('#'):
                # Extract sequence ID and length
                if 'Length:' in line:
                    parts = line.split()
                    seq_id = parts[1]
                    length = int(parts[3])
                    
                    current_sequence = {
                        'sequence_id': seq_id,
                        'length': length,
                        'tmr_count': 0,
                        'features': []
                    }
                    
                # Extract TMR count
                elif 'Number of predicted TMRs:' in line:
                    if current_sequence:
                        tmr_count = int(line.split(':')[-1].strip())
                        current_sequence['tmr_count'] = tmr_count
                        
            # Handle feature lines
            else:
                if current_sequence:
                    # Split by tab (GFF3 format)
                    fields = line.split('\t')
                    
                    if len(fields) >= 3:
                        feature = {
                            'seqid': fields[0],
                            'source': fields[1],
                            'type': fields[1],  # In this format, source and type are the same
                            'start': int(fields[2]),
                            'end': int(fields[3]) if len(fields) > 3 else int(fields[2]),
                            'score': fields[4] if len(fields) > 4 and fields[4] else None,
                            'strand': fields[5] if len(fields) > 5 and fields[5] else None,
                            'phase': fields[6] if len(fields) > 6 and fields[6] else None,
                            'attributes': fields[7] if len(fields) > 7 and fields[7] else None
                        }
                        current_sequence['features'].append(feature)
        
        # Don't forget the last sequence if file doesn't end with //
        if current_sequence:
            sequences.append(current_sequence)
    
    return sequences

def print_gff3_summary(sequences):
    """
    Print a summary of parsed GFF3 data.
    
    Args:
        sequences (list): Output from parse_gff3_tmr_file()
    """
    print(f"Parsed {len(sequences)} sequences:")
    print("-" * 50)
    
    for seq in sequences:
        print(f"Sequence: {seq['sequence_id']}")
        print(f"  Length: {seq['length']} aa")
        print(f"  Predicted TMRs: {seq['tmr_count']}")
        print(f"  Features: {len(seq['features'])}")
        
        for feature in seq['features']:
            print(f"    {feature['type']}: {feature['start']}-{feature['end']}")
        print()