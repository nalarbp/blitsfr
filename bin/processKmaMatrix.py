#!/usr/bin/env python
"""
Transform KMA matrix outputs into window-based coverage CGView plots.
It compiles per-base KMA matrix files into per-window coverage values and summary coverage tables.

Last checked by: BP
"""
import pandas as pd
import glob
import gzip
import os
import sys
import argparse
from typing import Dict, List, Tuple

def parse_args():
    """Parse arguments."""
    parser = argparse.ArgumentParser(description='Process KMA matrix files to create windowed coverage table')
    parser.add_argument('--kma_outputs', nargs='+', required=True, help='KMA matrix files (.mat.gz)')
    parser.add_argument('--res_file', required=True, help='Compiled KMA results file')
    parser.add_argument('--output', required=True, help='Output file for compiled matrix data')
    parser.add_argument('--coverage_output', required=True, help='Output file for coverage metrics (TSV format)')
    parser.add_argument('--output_normalised_1', help='Output file with values limited between min and max reads')
    parser.add_argument('--output_normalised_2', help='Output file with values normalised to range 0-1 (long format)')
    parser.add_argument('--window_size', type=int, default=1000, help='Window size in bp')
    parser.add_argument('--min_reads', type=float, default=0, help='Minimum read count (default: 0)')
    parser.add_argument('--max_reads', type=float, default=100, help='Maximum read count (default: 100)')
    
    return parser.parse_args()

def generate_kma_coverage(res_file: str) -> pd.DataFrame:
    """Get total covered reference bases for each sample from the KMA results."""
    try:
        res_df = pd.read_csv(res_file, sep='\t')
        
        #calculate bp_covered for each sample-template pair
        res_df['bp_covered'] = (res_df['Template_Coverage'] / 100) * res_df['Template_length']
        res_df['bp_covered'] = res_df['bp_covered'].round().astype(int)
        
        #group by SampleID and sum bp_covered
        coverage_results = res_df.groupby('SampleID')['bp_covered'].sum().reset_index()
        coverage_results = coverage_results.rename(columns={
            'SampleID': 'id',
            'bp_covered': 'total_bp_covered'
        })
        
        return coverage_results
        
    except Exception as e:
        print(f"Error generating KMA coverage: {e}", file=sys.stderr)
        return pd.DataFrame(columns=['query_file', 'total_bp_covered'])
    
    
def get_template_lengths_from_res(res_file: str) -> Dict[str, int]:
    """Collect template (ref.) lengths from the compiled KMA results."""
    template_lengths = {}
    try:
        res_df = pd.read_csv(res_file, sep='\t')
        for _, row in res_df.iterrows():
            template_id = row['#Template']
            length = row['Template_length']
            if template_id in template_lengths and template_lengths[template_id] != length:
                print(f"Warning: Template {template_id} has inconsistent lengths in res file", file=sys.stderr)
            template_lengths[template_id] = length
    except Exception as e:
        print(f"Error reading res file: {e}", file=sys.stderr)
    
    return template_lengths

def generate_windows(template_lengths: Dict[str, int], window_size: int) -> pd.DataFrame:
    """Generate fixed-width windows across each template sequence."""
    result_data = []
    
    for template_id, length in template_lengths.items():
        for start in range(1, length + 1, window_size):
            end = min(start + window_size - 1, length)
            window_name = f"window_{start}_{end}"
            result_data.append({
                'Template': template_id,
                'Window': window_name,
                'Start': start,
                'End': end
            })
    
    return pd.DataFrame(result_data)

def process_matrix_file(mat_file: str, window_size: int) -> Tuple[str, Dict[Tuple[str, str], float]]:
    """Convert one KMA matrix file into average coverage values per defined template window."""
    sample_id = os.path.basename(mat_file).replace('.mat.gz', '')
    sample_coverage = {}
    
    with gzip.open(mat_file, 'rt') as f:
        lines = f.readlines()
        
        current_template = None
        coverages = []
        position = 0
        
        for line_num, line in enumerate(lines):
            line = line.strip()
            if not line:
                continue
                
            #check header line
            if line.startswith('#'):
                if current_template and coverages:
                    #get average coverage for each window
                    for i in range(0, len(coverages), window_size):
                        start_pos = i + 1  #position is 1-based
                        end_pos = min(start_pos + window_size - 1, len(coverages))
                        window_coverage = coverages[i:min(i+window_size, len(coverages))]
                        if window_coverage:
                            avg_coverage = sum(window_coverage) / len(window_coverage)
                            window_name = f"window_{start_pos}_{end_pos}"
                            sample_coverage[(current_template, window_name)] = avg_coverage
                
                #when start on another new template
                current_template = line.replace('#', '')
                coverages = []
                position = 0
                continue
            
            parts = line.split()
            
            #check to have 7 columns
            if len(parts) < 7:
                print(f"Warning: Invalid data line in {mat_file} at line {line_num+1}: {line}", file=sys.stderr)
                continue
            
            #skip deletion positions
            if parts[0] == '-': #bug potential, what if actual ref has - character
                continue
            
            #calculate coverage (sum of A,C,G,T,N counts)
            try:
                coverage = sum(int(parts[i]) for i in range(1, 6))
                coverages.append(coverage)
            except (ValueError, IndexError) as e:
                print(f"Error processing line in {mat_file} at line {line_num+1}: {line} - {e}", file=sys.stderr)
                continue
        
        #on the last template
        if current_template and coverages:
            for i in range(0, len(coverages), window_size):
                start_pos = i + 1
                end_pos = min(start_pos + window_size - 1, len(coverages))
                window_coverage = coverages[i:min(i+window_size, len(coverages))]
                if window_coverage:
                    avg_coverage = sum(window_coverage) / len(window_coverage)
                    window_name = f"window_{start_pos}_{end_pos}"
                    sample_coverage[(current_template, window_name)] = avg_coverage
    
    return sample_id, sample_coverage

def normalise_dataframe(df: pd.DataFrame, min_reads: float, max_reads: float) -> pd.DataFrame:
    """Normalisation 1: clip per-window sample coverage values to the defined min/max reads range."""
    normalised_df = df.copy()
    
    #get the sample columns (no need Template, Window, Start, End)
    sample_columns = [col for col in df.columns if col not in ['Template', 'Window', 'Start', 'End']]
    
    for col in sample_columns:
        normalised_df[col] = normalised_df[col].clip(min_reads, max_reads)
    
    return normalised_df

def normalise_to_range(df: pd.DataFrame, min_reads: float, max_reads: float) -> pd.DataFrame:
    """Normalisation 2: scale windowed coverage values into a 0-1 long-format table for cgview.js plotting."""
    normalised_df = df.copy()
    sample_columns = [col for col in df.columns if col not in ['Template', 'Window', 'Start', 'End']]
    
    #check if using a single threshold (presence/absence mode)
    is_binary_mode = abs(max_reads - min_reads) < 1e-10  #using small epsilon for float comparison
    
    for col in sample_columns:
        if is_binary_mode:
            #binary presence/absence mode:
            #values >= threshold (min_reads/max_reads) become 1, values < threshold become 0
            threshold = min_reads
            normalised_df[col] = (normalised_df[col] >= threshold).astype(float)
        else:
            #standard continuous normalisation:
            normalised_df[col] = normalised_df[col].clip(min_reads, max_reads)
            #then normalise to 0-1 range
            normalised_df[col] = (normalised_df[col] - min_reads) / (max_reads - min_reads)
    
    #to long format
    long_df = pd.melt(
        normalised_df, 
        id_vars=['Template', 'Window', 'Start', 'End'],
        var_name='SampleID',
        value_name='Value'
    )
    #reorder columns
    long_df = long_df[['SampleID', 'Template', 'Window', 'Start', 'End', 'Value']]
    return long_df

def main():
    """Parse inputs, compile matrix windows, and write KMA output tables."""
    args = parse_args()
    window_size = args.window_size
    min_reads = args.min_reads
    max_reads = args.max_reads
    mat_files = args.kma_outputs if isinstance(args.kma_outputs, list) else glob.glob(args.kma_outputs)
    
    if abs(max_reads - min_reads) < 1e-10:
        print(f"Using binary presence/absence mode with threshold {min_reads}")
    
    template_lengths = get_template_lengths_from_res(args.res_file)
    master_df = generate_windows(template_lengths, window_size)
    
    #process each matrix file
    for mat_file in mat_files:
        sample_id, sample_coverage = process_matrix_file(mat_file, window_size)
        
        for idx, row in master_df.iterrows():
            key = (row['Template'], row['Window'])
            if key in sample_coverage:
                master_df.loc[idx, sample_id] = sample_coverage[key]
            else:
                master_df.loc[idx, sample_id] = 0
    
    #write out coverage
    if args.coverage_output:
        coverage_df = generate_kma_coverage(args.res_file)
        coverage_df.to_csv(args.coverage_output, sep='\t', index=False)
        print(f"Coverage metrics written to {args.coverage_output}")
    
    #write out compiled matrix data to a file
    master_df.to_csv(args.output, sep='\t', index=False)
    print(f"Successfully processed {len(mat_files)} matrix files. Output written to {args.output}")
    
    #normalised 1 out
    if args.output_normalised_1:
        normalised_df1 = normalise_dataframe(master_df, min_reads, max_reads)
        normalised_df1.to_csv(args.output_normalised_1, sep='\t', index=False)
        print(f"Normalised output (min/max limited) written to {args.output_normalised_1}")
    
    #normalised 2 out
    if args.output_normalised_2:
        normalised_df2 = normalise_to_range(master_df, min_reads, max_reads)
        
        if abs(max_reads - min_reads) < 1e-10:
            mode_desc = f"binary presence/absence (threshold: {min_reads})"
        else:
            mode_desc = f"continuous 0-1 scale (min: {min_reads}, max: {max_reads})"
        
        normalised_df2.to_csv(args.output_normalised_2, sep='\t', index=False)
        print(f"Normalised output (long format, {mode_desc}) written to {args.output_normalised_2}")

if __name__ == "__main__":
    main()
