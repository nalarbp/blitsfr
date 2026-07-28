#!/usr/bin/env python
"""
Calculate subject (ref,) coverage for each query track.
E.g. How many subject positions are covered by at least one alignment.

Last checked by: BP
"""
import argparse
import pandas as pd
import numpy as np
from concurrent.futures import ProcessPoolExecutor
import os

def process_group(group_data):
    """Calculate covered subject positions for each query-subject BLAST group."""
    query, subject = group_data[0]
    group = group_data[1]
    
    subject_length = group['slen'].iloc[0]
    coverage = np.zeros(subject_length, dtype=np.bool_)
    
    #set positions
    for _, row in group.iterrows():
        start = min(row['sstart'], row['send']) - 1
        end = max(row['sstart'], row['send'])
        coverage[start:end] = True
    
    #get total covered positions
    subject_positions_covered = np.sum(coverage)
    coverage_percentage = (subject_positions_covered / subject_length) * 100
    
    return {
        'query_file': query,
        'subject_id': subject,
        'subject_length': subject_length,
        'subject_positions_covered': subject_positions_covered,
        'subject_positions_covered_percentage': coverage_percentage
    }

def calculate_coverage(blast_df, num_cores=1):
    """Calculate coverage for all query-subject groups and return the final table."""
    #group by query_file and subject id (ref can have multiple subjects/contigs)
    groups = list(blast_df.groupby(['query_file', 'sseqid']))
    
    if num_cores > 1:
        try:
            with ProcessPoolExecutor(max_workers=num_cores) as executor:
                results = list(executor.map(process_group, groups))
        except PermissionError:
            # fall back to single-process mode when process pools are blocked
            results = [process_group(group_data) for group_data in groups]
    else:
        results = [process_group(group_data) for group_data in groups]
    
    results_df = pd.DataFrame(results)
    return results_df

def get_optimal_cores(threads_arg):
    """Just to check if the requested threads fit the machine capacity."""
    available_cores = os.cpu_count() or 1
    
    if threads_arg == 'max':
        return max(1, available_cores - 2)
    else:
        try:
            return min(int(threads_arg), available_cores)
        except ValueError:
            return max(1, available_cores - 2)

def main():
    """Parse args and write the BLAST coverage summary table."""
    parser = argparse.ArgumentParser(description='Process BLAST results to calculate coverage metrics')
    parser.add_argument('input_file', help='Input BLAST results file (TSV format)')
    parser.add_argument('output_file', help='Output file for coverage metrics (TSV format)')
    parser.add_argument('-t', '--threads', default='max',
                        help='Number of CPU threads to use. Can be an integer or "max" (default: max)')
    args = parser.parse_args()
    
    num_cores = get_optimal_cores(args.threads)
    
    try:
        blast_df = pd.read_csv(args.input_file, sep='\t')
        required_cols = ['query_file', 'sseqid', 'sstart', 'send', 'slen']
        missing_cols = [col for col in required_cols if col not in blast_df.columns]
        
        if missing_cols:
            raise ValueError(f"Missing required columns: {', '.join(missing_cols)}")
        
        print(f"Processing with {num_cores} CPU cores...")
        coverage_df = calculate_coverage(blast_df, num_cores)
        coverage_df = coverage_df.sort_values('subject_positions_covered', ascending=False)
        coverage_df.to_csv(args.output_file, sep='\t', index=False)
        print(f"Coverage metrics written to {args.output_file}")
    
    except Exception as e:
        print(f"Error: {e}")
        raise

if __name__ == "__main__":
    main()
