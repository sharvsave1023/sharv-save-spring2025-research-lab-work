#!/usr/bin/env python3
import time
from Bio import SeqIO
import pandas as pd
from riot_na import Riot

def load_fasta(fasta_file):
    """Load FASTA file and return sequences with IDs"""
    sequences = []
    ids = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append(str(record.seq))
        ids.append(record.id)
    return sequences, ids

def process_with_riot(fasta_file):
    print("Starting RIOT processing...")
    start_time = time.time()
    
    # Initialize RIOT
    riot = Riot()
    results = riot.annotate_sequences(fasta_file)
    
    # Convert RIOT results to DataFrame format
    processed_data = []
    for result in results:
        processed_data.append({
            'sequence_id': result['name'],
            'sequence': result['sequence'],
            'fr1': result['FR1'],
            'cdr1': result['CDR1'],
            'fr2': result['FR2'],
            'cdr2': result['CDR2'],
            'fr3': result['FR3'],
            'cdr3': result['CDR3'],
            'fr4': result['FR4'],
            'v_gene': result['V_gene'],
            'success': result['CDR3'] is not None
        })
    
    end_time = time.time()
    processing_time = end_time - start_time
    
    # Create DataFrame and calculate statistics
    df = pd.DataFrame(processed_data)
    success_rate = (df['success'].sum() / len(df)) * 100
    
    print(f"\nRIOT Processing Results:")
    print(f"Total sequences processed: {len(df)}")
    print(f"Processing time: {processing_time:.2f} seconds")
    print(f"Success rate: {success_rate:.2f}%")
    
    return df

def main():
    input_file = 'ABSD_1000.fasta'
    output_file = 'riot_results.csv'
    
    # Process sequences
    results_df = process_with_riot(input_file)
    
    # Save results
    results_df.to_csv(output_file, index=False)
    print(f"\nResults saved to {output_file}")

if __name__ == "__main__":
    main()