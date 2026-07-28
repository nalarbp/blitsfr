#!/usr/bin/env python
"""
Generate the final BLITSFR payload and inject it into the template.

Combines CGView JSON, metadata if given, coverage summaries, and pipeline log
into one payload that can be embedded into the final template report.

Last checked by: BP
"""
import argparse
import datetime
from pathlib import Path
import pandas as pd
import orjson
from scifrMutator import mutate_template_memory


def trim_display_path(value):
    """Trim absolute paths into shorter display paths."""
    trimmed = value.strip()
    if "/" not in trimmed and "\\" not in trimmed:
        return value
    return f"../{Path(trimmed).name}"


def sanitize_log_value(value):
    """Recursively trim path-like values."""
    if isinstance(value, dict):
        return {key: sanitize_log_value(val) for key, val in value.items()}
    if isinstance(value, list):
        return [sanitize_log_value(item) for item in value]
    if isinstance(value, str):
        if "," in value:
            parts = [part.strip() for part in value.split(",")]
            if all("/" in part or "\\" in part for part in parts if part):
                return ",".join(trim_display_path(part) for part in parts)
        return trim_display_path(value)
    return value

def combine_data(cgview_json, metadata_tsv, coverage_tsv, log_json, pipeline_version):
    """Build the embedded BLITSFR payload object."""
    #cgview
    cgview_data = 'NA'
    if cgview_json:
        try:
            with open(cgview_json, 'r') as cgview_json_file:
                cgview_data = orjson.loads(cgview_json_file.read())
        except Exception as e:
            print(f"Warning: Could not read cgview file: {str(e)}")
    
    #metadata
    metadata_data = 'NA'
    if metadata_tsv and metadata_tsv != 'false':
        try:
            with open(metadata_tsv, 'r') as metadata_tsv_file:
                tsv_content = metadata_tsv_file.read()
                metadata_data = tsv_content.replace("\t", ";t").replace("\n", ";n")
        except Exception as e:
            print(f"Warning: Could not read metadata file: {str(e)}")
    
    #coverage data
    coverage_data = 'NA'
    if coverage_tsv:
        try:
            coverage_df = pd.read_csv(coverage_tsv, sep='\t')
            
            #group by query_file (track) and sum the subject_positions_covered
            if 'query_file' in coverage_df.columns and 'subject_positions_covered' in coverage_df.columns:
                aggregated_coverage = coverage_df.groupby('query_file')['subject_positions_covered'].sum().reset_index()
                aggregated_coverage.columns = ['id', 'total_bp_covered']
                
                coverage_content = aggregated_coverage.to_csv(sep='\t', index=False)
                coverage_data = coverage_content.replace("\t", ";t").replace("\n", ";n")
            else:
                print("Warning: Coverage file doesn't have required columns 'query_file' and 'subject_positions_covered'")
                coverage_content = pd.read_csv(coverage_tsv, sep='\t').to_csv(sep='\t', index=False)
                coverage_data = coverage_content.replace("\t", ";t").replace("\n", ";n")
                
        except Exception as e:
            print(f"Warning: Could not process coverage file: {str(e)}")
    
    #log 
    log_data = {}
    if log_json:
        try:
            with open(log_json, 'r') as log_file:
                log_data = orjson.loads(log_file.read())
                log_data = sanitize_log_value(log_data)
                
        except Exception as e:
            print(f"Warning: Could not read log file: {str(e)}")
    
    log_data['pipeline_version'] = pipeline_version
    log_data["created"] = datetime.datetime.now().isoformat()
    
    combined_data = {
            "startIdx": "@@BLITSFR@@INPUT@@START@@",
            "logData": log_data,
            "trackMetadata": metadata_data,
            "coverageData": coverage_data,
            "cgviewData": cgview_data,
            "endIdx": "@@BLITSFR@@INPUT@@END@@"
    }
    
    return combined_data

def generate_scifr(cgview_json, metadata_json, coverage_json, template_path, output_path, json_output_path, log_json, pipeline_version, save_intermediate):
    """Build the final payload and write the mutated single-file report."""
    try:
        # Step 1: Combine data
        combined_data = combine_data(cgview_json, metadata_json, coverage_json, log_json, pipeline_version)
        
        # Step 2: Write JSON output
        if save_intermediate:
            json_bytes = orjson.dumps(combined_data, option=orjson.OPT_NON_STR_KEYS)
            with open(json_output_path, 'wb', buffering=1024*1024) as outfile:
                outfile.write(json_bytes)
        
        # Step 3: Process template
        mutate_template_memory(combined_data, template_path, output_path)
       
    except Exception as e:
        print(f"Error generating SCIFR: {str(e)}")
        raise

if __name__ == "__main__":
    """Parse CLI inputs and generate the final BLITSFR report."""
    parser = argparse.ArgumentParser(description="Generate SCIFR from JSONs")
    parser.add_argument("-c", "--cgview_json", help="CGview", required=True)
    parser.add_argument("-m", "--metadata", help="Metadata")
    parser.add_argument("-d", "--coverage", help="Coverage data")
    parser.add_argument("-t", "--template", help="Template file", required=True)
    parser.add_argument("-o", "--output", help="Output file", required=True)
    parser.add_argument("-j", "--json", help="Output JSON file", required=True)
    parser.add_argument("-l", "--log", help="Params JSON file", required=True)
    parser.add_argument("--save_intermediate", help="Save intermediate JSON", action='store_true', default=False)
    parser.add_argument("-p", "--pipeline_version", help="Pipeline version", required=True)
    
    args = parser.parse_args()
    generate_scifr(
        args.cgview_json,
        args.metadata,
        args.coverage,
        args.template,
        args.output,
        args.json,
        args.log,
        args.pipeline_version,
        args.save_intermediate
    )
