#!/usr/bin/env python3
import time
from Bio import SeqIO
import pandas as pd
from abnumber import Chain
from riot_na import Riot
from pandarallel import pandarallel

def load_fasta(fasta_file):
    """seq --> fasta --> dict"""
    sequences = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append(str(record.seq))
    return sequences

def process_with_abnumber(sequences):
    """seq thru abmumber"""
    results = []
    start_time = time.time()
    success_count = 0
    
    for seq in sequences:
        try:
            chain = Chain(seq, scheme='imgt', allowed_species='human', assign_germline=True)
            results.append({
                'sequence': seq,
                'fr1': chain.fr1_seq,
                'cdr1': chain.cdr1_seq,
                'fr2': chain.fr2_seq,
                'cdr2': chain.cdr2_seq,
                'fr3': chain.fr3_seq,
                'cdr3': chain.cdr3_seq,
                'v_gene': chain.v_gene,
                'success': True
            })
            success_count += 1
        except Exception as e:
            results.append({
                'sequence': seq,
                'success': False,
                'error': str(e)
            })
    
    end_time = time.time()
    processing_time = end_time - start_time
    
    return {
        'results': results,
        'time': processing_time,
        'success_rate': success_count / len(sequences)
    }

def process_with_riot(sequences):
    """seq thru RIOT"""
    start_time = time.time()
    
    # Create temporary FASTA file for RIOT
    with open('temp.fasta', 'w') as f:
        for i, seq in enumerate(sequences):
            f.write(f'>seq_{i}\n{seq}\n')
    
    # Initialize RIOT
    riot = Riot()
    results = riot.annotate_sequences('temp.fasta')
    
    end_time = time.time()
    processing_time = end_time - start_time
    
    # Count successful annotations
    success_count = sum(1 for result in results if result['CDR3'] is not None)
    
    return {
        'results': results,
        'time': processing_time,
        'success_rate': success_count / len(sequences)
    }

def main():
    # Initialize parallel processing for larger datasets
    pandarallel.initialize()
    
    # Load sequences from FASTA file
    print("Loading sequences from FASTA file...")
    sequences = load_fasta('absd.fasta')
    total_sequences = len(sequences)
    print(f"Loaded {total_sequences} sequences")
    
    # Process with Abnumber
    print("\nProcessing with Abnumber...")
    abnumber_results = process_with_abnumber(sequences)
    
    # Process with RIOT
    print("\nProcessing with RIOT...")
    riot_results = process_with_riot(sequences)
    
    # Print comparison results
    print("\nPerformance Comparison:")
    print("-" * 50)
    print(f"Total sequences processed: {total_sequences}")
    print("\nAbnumber:")
    print(f"Processing time: {abnumber_results['time']:.2f} seconds")
    print(f"Success rate: {abnumber_results['success_rate']*100:.2f}%")
    
    print("\nRIOT:")
    print(f"Processing time: {riot_results['time']:.2f} seconds")
    print(f"Success rate: {riot_results['success_rate']*100:.2f}%")
    
    # Save detailed results to CSV files
    pd.DataFrame(abnumber_results['results']).to_csv('abnumber_results.csv', index=False)
    pd.DataFrame(riot_results['results']).to_csv('riot_results.csv', index=False)

if __name__ == "__main__":
    main()