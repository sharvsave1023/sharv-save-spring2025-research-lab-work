#!/usr/bin/env python3
import pandas as pd
import numpy as np
from Bio import pairwise2
from pandarallel import pandarallel
import time
from Levenshtein import distance as lev_distance

def load_results():
    """Load results from both RIOT and ABNumber CSV files"""
    print("Loading results...")
    riot_df = pd.read_csv('riot_results.csv')
    abn_df = pd.read_csv('abnumber_results.csv')
    
    riot_df = riot_df.sort_values('sequence_id').reset_index(drop=True)
    abn_df = abn_df.sort_values('sequence_id').reset_index(drop=True)
    
    return riot_df, abn_df

def calculate_region_differences(row, region):
    """Calculate Levenshtein distance between RIOT and ABNumber regions"""
    riot_seq = str(row[f'riot_{region}']) if pd.notna(row[f'riot_{region}']) else ''
    abn_seq = str(row[f'abn_{region}']) if pd.notna(row[f'abn_{region}']) else ''
    
    if riot_seq == '' and abn_seq == '':
        return np.nan
    elif riot_seq == '' or abn_seq == '':
        return -1 
    else:
        return lev_distance(riot_seq, abn_seq)

def compare_annotations():
    pandarallel.initialize(progress_bar=True)
    riot_df, abn_df = load_results()
    
    regions = ['fr1', 'cdr1', 'fr2', 'cdr2', 'fr3', 'cdr3', 'fr4']
    
    comparison_df = pd.DataFrame()
    comparison_df['sequence_id'] = riot_df['sequence_id']
    
    for region in regions:
        comparison_df[f'riot_{region}'] = riot_df[region]
        comparison_df[f'abn_{region}'] = abn_df[region]
    
    # Calculate Levenshtein distances in parallel
    print("\nCalculating Levenshtein distances...")
    for region in regions:
        comparison_df[f'{region}_distance'] = comparison_df.parallel_apply(
            lambda row: calculate_region_differences(row, region), axis=1
        )
    
    return comparison_df

def analyze_differences(comparison_df):
    """Analyze and display differences between RIOT and ABNumber results"""
    regions = ['fr1', 'cdr1', 'fr2', 'cdr2', 'fr3', 'cdr3', 'fr4']
    
    print("\nAnalysis of differences:")
    print("-" * 50)
    
    # For each region, analyze differences
    for region in regions:
        distances = comparison_df[f'{region}_distance']
        total_comparisons = len(distances.dropna())
        different = len(distances[distances > 0])
        failed_comparisons = len(distances[distances == -1])
        identical = len(distances[distances == 0])
        
        print(f"\n{region.upper()} Analysis:")
        print(f"Total comparisons: {total_comparisons}")
        print(f"Identical annotations: {identical}")
        print(f"Different annotations: {different}")
        print(f"Failed comparisons: {failed_comparisons}")
        
        if different > 0:
            avg_distance = distances[distances > 0].mean()
            print(f"Average Levenshtein distance (when different): {avg_distance:.2f}")
    
    # Save detailed differences to CSV
    save_detailed_differences(comparison_df)

def save_detailed_differences(comparison_df):
    """Save detailed differences to CSV file"""
    regions = ['fr1', 'cdr1', 'fr2', 'cdr2', 'fr3', 'cdr3', 'fr4']
    
    # Create mask for rows with any differences
    diff_mask = False
    for region in regions:
        diff_mask |= (comparison_df[f'{region}_distance'] > 0)
    
    # Filter and save differences
    differences_df = comparison_df[diff_mask]
    differences_df.to_csv('annotation_differences.csv', index=False)
    print(f"\nDetailed differences saved to annotation_differences.csv")
    
    # Create summary table of differences
    summary_rows = []
    for idx, row in differences_df.iterrows():
        for region in regions:
            if row[f'{region}_distance'] > 0:
                summary_rows.append({
                    'sequence_id': row['sequence_id'],
                    'region': region,
                    'RIOT': row[f'riot_{region}'],
                    'ABNumber': row[f'abn_{region}'],
                    'Levenshtein_distance': row[f'{region}_distance']
                })
    
    pd.DataFrame(summary_rows).to_csv('difference_summary.csv', index=False)
    print(f"Difference summary saved to difference_summary.csv")

def main():
    start_time = time.time()
    
    # Compare annotations
    comparison_df = compare_annotations()
    
    # Analyze and display differences
    analyze_differences(comparison_df)
    
    end_time = time.time()
    print(f"\nTotal processing time: {end_time - start_time:.2f} seconds")

if __name__ == "__main__":
    main()