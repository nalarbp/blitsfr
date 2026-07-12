#!/usr/bin/env python
import argparse
import os
import subprocess
import sys
import glob
from pathlib import Path

def transform_inputs_for_nextflow(inputs):
    #multiple input files, comma-separated 
    if ',' in inputs:
        files = [f.strip().strip('"\'') for f in inputs.split(',')]
        expanded_files = []
        for f in files:
            if '*' in f or '?' in f or '[' in f:
                matches = glob.glob(f)
                expanded_files.extend(matches)
            else:
                expanded_files.append(f)
        
        #return comma-separated
        return ','.join(expanded_files)
    #else as is (single or glob)
    else:
        return inputs
    
def run_blitsfr_pipeline(version, reference, method, queries, metadata, output, 
                         title, cpu_per_task,
                         precluster_queries, precluster_queries_args,
                         resume, config, nf_args, **method_params):
    if not os.path.exists(reference):
        sys.exit(f"Error: Reference file {reference} does not exist.")

    #determine method-specific parameters
    if method == 'blast':
        queries_fasta = transform_inputs_for_nextflow(queries)
        queries_fastq = "false"
    else:  # method == 'kma'
        queries_fasta = "false"
        queries_fastq = queries 
        
        #validate KMA reads mode
        kma_reads_mode = method_params.get('reads_mode', 'paired')
        if kma_reads_mode not in ['single', 'paired']:
            sys.exit(f"Error: KMA reads mode must be either 'single' or 'paired', got {kma_reads_mode}")
        
        #handle single-end reads
        if kma_reads_mode == 'single':
            queries_fastq = transform_inputs_for_nextflow(queries_fastq)
            
        #handle missing paired
        if kma_reads_mode == 'paired' and '{1,2}' not in queries_fastq and 'R{1,2}' not in queries_fastq:
                sys.exit(f"Error: Paired-end mode requires paired files. Use pattern '*_R{{1,2}}*.gz' to capture both files.")
                #handle wildcard in queries
                if '*' in queries_fastq:
                    if kma_reads_mode == 'paired' and '{1,2}' not in queries_fastq and 'R{1,2}' not in queries_fastq:
                        sys.exit(f"Error: Paired-end mode requires a pattern like '*_R{{1,2}}*.gz' to capture both files.")
                        
                    query_files = glob.glob(queries_fastq.replace('{1,2}', '1'))  # Just check for _1 files as test
                    if not query_files and '{1,2}' in queries_fastq:
                        sys.exit(f"Error: No read files found matching {queries_fastq}")
                elif not os.path.exists(queries_fastq):
                    sys.exit(f"Error: Reads path {queries_fastq} does not exist.")

    if metadata != "false" and not os.path.exists(metadata):
        sys.exit(f"Error: Metadata file {metadata} does not exist.")

    if precluster_queries:
        print(
            "Warning: --precluster_queries is under development and is currently disabled internally. Proceeding with preclustering turned off.",
            file=sys.stderr,
        )
        precluster_queries = False
        precluster_queries_args = ""

    #find the main.nf script
    script_path = os.path.realpath(os.path.abspath(__file__))
    script_dir = os.path.dirname(script_path)
    conda_prefix = os.environ.get('CONDA_PREFIX')

    if conda_prefix and os.path.exists(os.path.join(conda_prefix, 'blitsfr', 'main.nf')):
        script_dir = os.path.join(conda_prefix, 'blitsfr')

    main_nf_path = os.path.join(script_dir, 'main.nf')
    if not os.path.exists(main_nf_path):
        parent_dir = os.path.dirname(script_dir)
        main_nf_path = os.path.join(parent_dir, 'main.nf')

    if not os.path.exists(main_nf_path):
        sys.exit("Error: Cannot find main.nf. Please ensure BLITSFR is correctly installed.")

    #build the nextflow command
    print(f"Using workflow: {main_nf_path}")
    nextflow_cmd = ["nextflow", "run", main_nf_path]
    
    #add version
    nextflow_cmd.append(f"--pipeline_version={version}")

    #add method parameter
    nextflow_cmd.append(f"--method={method}")
    
    #add reference
    nextflow_cmd.append(f"--reference_genbank={reference}")
    
    #add queries based on method
    if method == 'blast':
        nextflow_cmd.append(f"--queries_fasta={queries_fasta}")
        nextflow_cmd.append("--queries_fastq=false")
        nextflow_cmd.append(f"--blast_filter_min_identity={method_params.get('blast_filter_min_identity', 80)}")
        nextflow_cmd.append(f"--blast_filter_min_coverage={method_params.get('blast_filter_min_coverage', 0)}")
        nextflow_cmd.append(f"--blast_filter_min_alignment={method_params.get('blast_filter_min_alignment', 200)}")
        nextflow_cmd.append(f"--precluster_queries={'true' if precluster_queries else 'false'}")
        
        #add BLAST-specific parameters
        blast_args = method_params.get('blast_args', '')
        if blast_args:
            nextflow_cmd.append(f"--blast_args={blast_args}")
        if precluster_queries_args:
            nextflow_cmd.append(f"--precluster_queries_args='{precluster_queries_args}'")
    else:  #method == 'kma'
        nextflow_cmd.append(f"--queries_fastq={queries_fastq}")
        nextflow_cmd.append("--queries_fasta=false")
        
        #add KMA-specific parameters
        nextflow_cmd.append(f"--kma_window={method_params.get('window', 1000)}")
        nextflow_cmd.append(f"--kma_reads_mode={method_params.get('reads_mode', 'paired')}")
        nextflow_cmd.append(f"--kma_min_reads={method_params.get('min_reads', 0)}")
        nextflow_cmd.append(f"--kma_max_reads={method_params.get('max_reads', 100)}")
        
        kma_args = method_params.get('kma_args', '')
        if kma_args:
            nextflow_cmd.append(f"--kma_args='{kma_args}'")

    if metadata != "false":
        nextflow_cmd.append(f"--tracks_metadata={metadata}")

    nextflow_cmd.extend([
        f"--results_directory={output}",
        f"--cgview_title='{title}'",
        f"--cpu_per_task={cpu_per_task}"
    ])

    if resume:
        nextflow_cmd.append("-resume")

    if config:
        nextflow_cmd.append(f"-c {config}")

    if nf_args:
        nextflow_cmd.append(nf_args)

    #run the command
    print(f"Executing: {' '.join(nextflow_cmd)}")
    try:
        subprocess.run(nextflow_cmd, check=True)
        print("\nBLITSFR completed successfully!")
    except subprocess.CalledProcessError as e:
        sys.exit(f"\nBLITSFR failed with exit code {e.returncode}")
    except KeyboardInterrupt:
        print("\nBLITSFR execution interrupted by user")
        sys.exit(1)

def setup_common_arguments(parser):
    parser.add_argument("-r", "--reference", required=True, 
                        help="Reference sequence file in GenBank format")
    parser.add_argument("-m", "--metadata", default="false", 
                        help="Metadata file (TSV) to be associated to the each track (query). Must have column 'id' correspoding to the query file basename (Without file extension, e.g Sample_1.fna will have id: Sample_1)")
    parser.add_argument("-o", "--output", default="results", 
                        help="Output directory for results")
    parser.add_argument("--title", default='CGView Map', 
                        help="Title for the CGView visualization")
    parser.add_argument("--cpu_per_task", type=int, default=2, 
                        help="Number of CPU per nextflow TASK to use. Give this arg to 2 CPUs, on a system with total 8 CPUs will allow maximum of 3 parralel tasks. The 1-2 CPUs are reseved for headroom")
    
    # Nextflow-specific options
    parser.add_argument("--resume", action="store_true", 
                        help="Resume previous run (Nextflow -resume flag)")
    parser.add_argument("-c", "--config", default=None, 
                        help="Nextflow configuration file to use")
    parser.add_argument("--nf-args", default="", 
                        help="Additional arguments to pass to Nextflow (as a quoted string)")

def main():
    parser = argparse.ArgumentParser(
        prog="blitsfr",
        description="BLITSFR can be run on assemblies or reads mode. Try: blitsfr assemblies -h or blitsfr reads -h",
        epilog="BLITSFR: BLAST Interactive Tracks in a Single File Report. By Budi Permana (nalar.bp@gmail.com)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    #version
    version = '0.1.2'
    parser.add_argument("-v", "--version", action="version", version=f"%(prog)s {version}")
    
    #subparsers for modes
    subparsers = parser.add_subparsers(dest="command", help="Analysis type to perform")
    subparsers.required = True
    
    #assemblies (BLAST) subcommand
    blast_parser = subparsers.add_parser("assemblies", 
                                        help="Perform BLASTn between query DNA sequence (in fasta format) against a reference sequence (in genbank format)",
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                        epilog="Example: blitsfr assemblies -r reference.gbk -q 'assemblies/*.fasta' -o results")
    blast_parser.add_argument("-q", "--queries", required=True,
                            help="Path to FASTA file(s) for assembly queries. Can be a single file or multiple files using pattern like 'dir/*.fna' (IMPORTANT!: must be given as a quoted string)")
    blast_parser.add_argument("--blast-args", default="-dust no -evalue 1E-20", 
                            help="Additional arguments to pass to BLAST (as a quoted string)")
    blast_parser.add_argument("--blast-filter-min-identity", type=float, default=80, 
                        help="Minimum identity percentage for FILTERING step on BLAST results for CGView")
    blast_parser.add_argument("--blast-filter-min-coverage", type=float, default=0, 
                        help="Minimum coverage percentage (alignment_length / subject_length) * 100) for FILTERING step on BLAST results for CGView")
    blast_parser.add_argument("--blast-filter-min-alignment", type=float, default=200, 
                        help="Minimum alignment length (in bp) for FILTERING step on BLAST results for CGView")
    blast_parser.add_argument("--precluster_queries", action="store_true",
                            help="Under development. Intended to enable skani-based preclustering of assembly queries before report decomposition. Currently disabled internally and defaults to false.")
    blast_parser.add_argument("--precluster_queries_args", default="",
                            help="Under development. Additional arguments for query preclustering. Currently ignored because preclustering is disabled internally.")
    
    setup_common_arguments(blast_parser)
    
    #reads (KMA) subcommand
    kma_parser = subparsers.add_parser("reads", 
                                      help="Perform KMA read mapping between query reads (in fastq format) against a reference sequence (in genbank format)",
                                      formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                      epilog="Example: blitsfr reads --reads-mode paired -r reference.gbk -q 'reads/*_R{1,2}.fastq.gz' -o results")
    kma_parser.add_argument("-q", "--queries", required=True,
                          help="Path to FASTQ file(s) for read queries. For paired-end mode, use a pattern like 'dir/*_R{1,2}.fastq.gz' (IMPORTANT!: must be given as a quoted string)")
    kma_parser.add_argument("--reads-mode", choices=['single', 'paired'], default='paired',
                          help="Read mode for KMA (single or paired)")
    kma_parser.add_argument("--window", type=int, default=1000, 
                          help="Window size in bp for KMA score averaging. Window size of 1000, will average the number of reads covering each bp in the 1000-bp window.")
    kma_parser.add_argument("--min-reads", type=float, default=0, 
                          help="Minimum read count for CGView plot data normalisation")
    kma_parser.add_argument("--max-reads", type=float, default=100, 
                          help="Maximum read count for CGView plot data normalisation")
    kma_parser.add_argument("--kma-args", default="", 
                          help="Additional arguments to pass to KMA (as a quoted string)")
    setup_common_arguments(kma_parser)
    
    args = parser.parse_args()
    
    #common parameters for all modes
    common_params = {
        'reference': args.reference,
        'metadata': args.metadata,
        'output': args.output,
        'title': args.title,
        'cpu_per_task': args.cpu_per_task,
        'precluster_queries': getattr(args, 'precluster_queries', False),
        'precluster_queries_args': getattr(args, 'precluster_queries_args', ''),
        'resume': args.resume,
        'config': args.config,
        'nf_args': args.nf_args,
        'queries': args.queries
    }
    
    if args.command == "assemblies":
        method = "blast"
        method_params = {
            'blast_args': args.blast_args,
            'blast_filter_min_identity': args.blast_filter_min_identity,
            'blast_filter_min_coverage': args.blast_filter_min_coverage,
            'blast_filter_min_alignment': args.blast_filter_min_alignment,
        }
    elif args.command == "reads":
        method = "kma"
        method_params = {
            'reads_mode': args.reads_mode,
            'window': args.window,
            'kma_args': args.kma_args,
            'min_reads': args.min_reads,
            'max_reads': args.max_reads,
        }
    else:
        sys.exit(f"Unknown command: {args.command}")
    
    #run the nextflow pipeline with the appropriate parameters
    run_blitsfr_pipeline(
        version=version,
        method=method,
        **common_params,
        **method_params
    )

if __name__ == "__main__":
    main()
