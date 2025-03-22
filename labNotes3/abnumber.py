#!/usr/bin/env python3
import time
from Bio import SeqIO
import pandas as pd
from abnumber import Chain

def load_fasta(fasta_file):
    """Load FASTA file and return sequences with IDs"""
    sequences = []
    ids = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append(str(record.seq))
        ids.append(record.id)
    return sequences, ids

def process_with_abnumber(sequences, ids):
    print("Starting Abnumber processing...")
    start_time = time.time()
    
    processed_data = []
    for seq, seq_id in zip(sequences, ids):
        try:
            chain = Chain(seq, scheme='imgt', allowed_species='human', assign_germline=True)
            processed_data.append({
                'sequence_id': seq_id,
                'sequence': seq,
                'fr1': chain.fr1_seq,
                'cdr1': chain.cdr1_seq,
                'fr2': chain.fr2_seq,
                'cdr2': chain.cdr2_seq,
                'fr3': chain.fr3_seq,
                'cdr3': chain.cdr3_seq,
                'fr4': chain.fr4_seq,
                'v_gene': chain.v_gene,
                'success': True
            })
        except Exception as e:
            processed_data.append({
                'sequence_id': seq_id,
                'sequence': seq,
                'fr1': None,
                'cdr1': None,
                'fr2': None,
                'cdr2': None,
                'fr3': None,
                'cdr3': None,
                'fr4': None,
                'v_gene': None,
                'success': False,
                'error': str(e)
            })
    
    end_time = time.time()
    processing_time = end_time - start_time
    
    # Create DataFrame and calculate statistics
    df = pd.DataFrame(processed_data)
    success_rate = (df['success'].sum() / len(df)) * 100
    
    print(f"\nAbnumber Processing Results:")
    print(f"Total sequences processed: {len(df)}")
    print(f"Processing time: {processing_time:.2f} seconds")
    print(f"Success rate: {success_rate:.2f}%")
    
    return df

def main():
    input_file = 'ABSD_1000.fasta'
    output_file = 'abnumber_results.csv'
    
    # Load sequences
    sequences, ids = load_fasta(input_file)
    
    # Process sequences
    results_df = process_with_abnumber(sequences, ids)
    
    # Save results
    results_df.to_csv(output_file, index=False)
    print(f"\nResults saved to {output_file}")

if __name__ == "__main__":
    main()