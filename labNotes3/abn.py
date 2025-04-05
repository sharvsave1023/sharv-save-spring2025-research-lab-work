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
    print("\nFirst few sequences:")
    for i in range(min(4, len(sequences))):
        print(f"Sequence {i}: ID={ids[i]}")
        print(f"Sequence: {sequences[i][:50]}...")
    
    # pair seq.
    for i in range(0, len(sequences), 2):
        try:
            hc_seq = sequences[i]
            lc_seq = sequences[i + 1] if (i + 1) < len(sequences) else None
            hc_id = ids[i]
            lc_id = ids[i + 1] if (i + 1) < len(ids) else None
            
            # Debug print
            if i < 4:
                print(f"\nProcessing pair {i//2}:")
                print(f"HC ID: {hc_id}")
                print(f"HC seq: {hc_seq[:50]}...")
                print(f"LC ID: {lc_id}")
                print(f"LC seq: {lc_seq[:50] if lc_seq else 'None'}...")
            
            # HC
            try:
                hc_chain = Chain(hc_seq, scheme='imgt', allowed_species='human', assign_germline=True)
            except Exception as e:
                print(f"\nError processing heavy chain {i//2}:")
                print(f"Sequence: {hc_seq[:50]}...")
                print(f"Error: {str(e)}")
                raise e
            
            hc_data = {
                'HC.Info': hc_id,
                'HC.seq': hc_seq,
                'LC.Info': lc_id,
                'LC.seq': lc_seq,
                'fr1': hc_chain.fr1_seq,
                'cdr1': hc_chain.cdr1_seq,
                'fr2': hc_chain.fr2_seq,
                'cdr2': hc_chain.cdr2_seq,
                'fr3': hc_chain.fr3_seq,
                'cdr3': hc_chain.cdr3_seq,
                'fr4': hc_chain.fr4_seq,
                'v_gene': hc_chain.v_gene,
                'success': True
            }
            
            # LC
            if lc_seq:
                try:
                    lc_chain = Chain(lc_seq, scheme='imgt', allowed_species='human', assign_germline=True)
                except Exception as e:
                    print(f"\nError processing light chain {i//2}:")
                    print(f"Sequence: {lc_seq[:50]}...")
                    print(f"Error: {str(e)}")
                    raise e
                
                hc_data.update({
                    'LC.Info': lc_id,
                    'LC.seq': lc_seq,
                    'lfr1': lc_chain.fr1_seq,
                    'lcdr1': lc_chain.cdr1_seq,
                    'lfr2': lc_chain.fr2_seq,
                    'lcdr2': lc_chain.cdr2_seq,
                    'lfr3': lc_chain.fr3_seq,
                    'lcdr3': lc_chain.cdr3_seq,
                    'lfr4': lc_chain.fr4_seq,
                    'light_v_gene': lc_chain.v_gene
                })
            
            processed_data.append(hc_data)
            
        except Exception as e:
            if i < 4: 
                print(f"\nError processing pair {i//2}:")
                print(f"Error: {str(e)}")
            processed_data.append({
                'HC.Info': hc_id,
                'HC.seq': hc_seq,
                'LC.Info': lc_id,
                'LC.seq': lc_seq,
                'success': False,
                'error': str(e)
            })
    
    end_time = time.time()
    processing_time = end_time - start_time
    
    df = pd.DataFrame(processed_data)
    success_rate = (df['success'].sum() / len(df)) * 100
    
    print(f"\nAbnumber Processing Results:")
    print(f"Total sequences processed: {len(df)}")
    print(f"Processing time: {processing_time:.2f} seconds")
    print(f"Success rate: {success_rate:.2f}%")
    
    if success_rate == 0 and 'error' in df.columns:
        print("\nSample of errors encountered:")
        print(df['error'].value_counts().head())
    return df

def main():
    input_file = 'absd_1000.fasta'
    output_file = 'abnumber_results.csv'
    
    sequences, ids = load_fasta(input_file)
    print(f"Loaded {len(sequences)} sequences")
    
    results_df = process_with_abnumber(sequences, ids)
    
    results_df.to_csv(output_file, index=False)
    print(f"\nResults saved to {output_file}")

if __name__ == "__main__":
    main()