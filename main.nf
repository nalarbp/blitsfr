#!/usr/bin/env nextflow

//Params
//Input output
params.reference_genbank = null
params.queries_fasta = null
params.queries_fastq = null
params.tracks_metadata = null
params.results_directory = null
params.method = null 
params.blast_filter_min_identity = null
params.blast_filter_min_coverage = null
params.blast_filter_min_alignment = null
params.blast_args = null
params.kma_window = null
params.kma_reads_mode = null
params.kma_args = ''
params.kma_min_reads = null
params.kma_max_reads = null
params.precluster_queries = false
params.precluster_queries_args = ''
params.max_parallel_jobs = Math.max(1, ((Runtime.getRuntime().availableProcessors() - 1) / params.cpu_per_task).intValue())
params.cgview_title = null
params.pipeline_version = null
params.save_intermediate_files = false
params.executor = "local"
params.cpu_per_task = 2

//Const
def BLAST_FORMAT = 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qcovs qcovhsp sstrand slen'
def MIN_REQUIRED_CPUS = params.cpu_per_task + 1
def PRECLUSTER_ENABLED = params.precluster_queries in [true, 'true']

//Resource validation
if (params.executor == "local" && Runtime.getRuntime().availableProcessors() < MIN_REQUIRED_CPUS) {
    log.error """
    ERROR: Insufficient CPU resources.

    This pipeline requires at least ${MIN_REQUIRED_CPUS} CPU cores with current settings:
    - ${params.cpu_per_task} cores for each task
    - 1 cores reserved

    Your system has ${Runtime.getRuntime().availableProcessors()} CPU cores available.

    Options:
    1. Run on a machine with more CPU cores
    2. Reduce CPU requirements by setting --cpu_per_task 1
    """
    exit 1
}

//Input validation
if (params.reference_genbank == null) {
    error("ERROR: Required parameter 'reference_genbank' not specified!")
}

if (params.method == 'blast' && params.queries_fasta == null) {
    error("ERROR: Required parameter 'queries_fasta' not specified for BLAST method!")
}

if (params.method == 'kma') {
    if (params.queries_fastq == null) {
        error("ERROR: Required parameter 'queries_fastq' not specified for KMA method!")
    }

    if (params.kma_reads_mode != "single" && params.kma_reads_mode != "paired") {
        error("ERROR: Invalid KMA reads mode! Must be either 'single' or 'paired'.")
    }
}

if (PRECLUSTER_ENABLED && params.method != 'blast') {
    error("ERROR: precluster_queries is currently only supported for BLAST assemblies mode!")
}

if (params.queries_fasta == null && params.queries_fastq == null) {
    error("ERROR: Either 'queries_fasta' or 'queries_fastq' must be specified!")
}

//Method validation
if (params.method != 'blast' && params.method != 'kma') {
    error("ERROR: Method must be either 'blast' or 'kma'!")
}

//Workflow
workflow {

    ch_ref = Channel.fromPath(params.reference_genbank, checkIfExists: true)
    ref_files = SETUP_REFERENCE(ch_ref)

    if (params.method == 'blast') {
        //BLAST workflow
        ref_db_files = MAKE_BLASTDB(ref_files.ref_fna).blast_db_files.collect()

        ch_query_files = Channel
            .fromPath(
                params.queries_fasta.contains(',') ? 
                params.queries_fasta.split(',').collect { it.trim().replaceAll('"', '') } : 
                params.queries_fasta.trim().replaceAll('"', ''), 
            checkIfExists: true)
            .map { file -> return tuple(file.baseName, file) }

        ch_query = ch_query_files.combine(ref_db_files)

        blast_results = RUN_BLAST(ch_query).blast_out
        all_blast_results = blast_results.collect()
        compiled_results = COMPILE_BLAST_RESULTS(all_blast_results)

        if (PRECLUSTER_ENABLED) {
            PRECLUSTER_QUERIES(
                ch_query_files.map { query_name, query_fasta -> query_fasta }.collect(),
                ref_files.ref_fna,
            )
        }

        ch_meta = params.tracks_metadata == "false"
            ? Channel.of([])
            : Channel.fromPath(params.tracks_metadata, checkIfExists: true)

        cgview_json = BUILD_CGVIEW_JSON_BLAST(
            compiled_results.compiled_results,
            ref_files.ref_fna,
            ref_files.ref_features,
        )

        ch_coverage = compiled_results.blast_coverage
    } else {
        //KMA workflow
        ref_kma_db_files = INDEX_KMADB(ref_files.ref_fna).kma_db_files.collect()
        if (params.kma_reads_mode == "paired") {
            log.info "Running in paired-end reads mode for KMA"
            ch_query_fastq = Channel
                .fromFilePairs("${params.queries_fastq}", checkIfExists: true)
                .ifEmpty {
                    error """
                    No paired fastq files found with the pattern: ${params.queries_fastq}

                    For paired-end mode, please use a pattern that captures read pairs such as:
                    --queries_fastq "path/to/reads/*_R{1,2}*"

                    Examples:
                    --queries_fastq "data/*_R{1,2}.fastq.gz"
                    --queries_fastq "data/*_R{1,2}.fq.gz"
                    """
                }

            //Combine with KMA database and run KMA
            ch_query_kma_input = ch_query_fastq
                .combine(ref_kma_db_files)
                
            kma_results = RUN_KMA_PAIRED(ch_query_kma_input)
        } else {
            log.info "Running in single-end reads mode for KMA"
            ch_query_fastq = Channel
                .fromPath("${params.queries_fastq}", checkIfExists: true)
                .map { file -> return tuple(file.baseName, file) }
                .combine(ref_kma_db_files)

            kma_results = RUN_KMA_SINGLE(ch_query_fastq)
        }

        compiled_results = COMPILE_KMA_RESULTS(kma_results.kma_res.collect(), 
            kma_results.kma_mat_gz.collect())

        ch_meta = params.tracks_metadata == "false"
            ? Channel.of([])
            : Channel.fromPath(params.tracks_metadata, checkIfExists: true)

        cgview_json = BUILD_CGVIEW_JSON_KMA(
            compiled_results.compiled_kma_mat_norm_2,
            ref_files.ref_fna,
            ref_files.ref_features,
        )

        ch_coverage = compiled_results.kma_coverage
    }

    template_file = file("${baseDir}/bin/blitsfr_template.html")
    scifr_output = GENERATE_SCIFR_REPORT(
        cgview_json,
        ch_meta,
        ch_coverage,
        template_file,
    )

    VALIDATE_REPORT(scifr_output.scifr_report)
}

//Analytical processes
process SETUP_REFERENCE {
    publishDir "${params.results_directory}/1_reference", mode: 'copy'

    input:
    path reference

    output:
    path "ref.fna", emit: ref_fna
    path "ref_features.gff3", emit: ref_features

    script:
    """
    # Setup reference files
    setupReference.py --ref_gbk ${reference} \
        --out_seqs ref.fna \
        --out_features ref_features.gff3
    """
}

process MAKE_BLASTDB {
    publishDir "${params.results_directory}/1_reference/blast_db", mode: 'copy'

    input:
    path ref_fna

    output:
    path "ref_db.*", emit: blast_db_files

    script:
    """
    # Create BLAST database from the reference fasta
    makeblastdb \
        -in ${ref_fna} \
        -dbtype nucl \
        -out ref_db
    """
}

process INDEX_KMADB {
    publishDir "${params.results_directory}/1_reference/kma_db", mode: 'copy'

    input:
    path ref_fna

    output:
    path "kma_db.*", emit: kma_db_files

    script:
    """
    # Create KMA database from the reference fasta
    kma index -i ${ref_fna} -o kma_db
    """
}

process RUN_BLAST {
    publishDir "${params.results_directory}/2_blast", mode: 'copy'

    input:
    tuple val(query_name), path(query_fasta), path("ref_db.ndb"), path("ref_db.nhr"), path("ref_db.nin"), path("ref_db.njs"), path("ref_db.not"), path("ref_db.nsq"), path("ref_db.ntf"), path("ref_db.nto")

    output:
    path "${query_name}.blast.out", emit: blast_out
    path "${query_name}.blast.out.raw", optional: true, emit: blast_raw

    script:
    def min_identity = params.blast_filter_min_identity ? params.blast_filter_min_identity : "80"
    def min_coverage = params.blast_filter_min_coverage ? params.blast_filter_min_coverage : "0"
    def min_alignment = params.blast_filter_min_alignment ? params.blast_filter_min_alignment : "100"
    def save_raw = params.save_intermediate_files ? '' : '&& rm ${query_name}.blast.out.raw'
    """
    # Run blast
    # Bug alert: subject_besthit or max_hsps = 1, will choose the the best hit.
    # What happend if contigs (query) is mapped to 2 different region in the subject
    # With those pramas set true, it will just choose 1 (the best)
    # Example is when looking at pUTI seq vs assembly when max_hsps =1 . somehow the replicon identitdied 
    # by plasmidfinder has no hits displayed, its because the hit on replicon is from the node (contig)
    # that has another hits, which better hits (longer aln, etc) to other region.
    # Not setting max_hsps is gurantee to see all hits (but is this biologically relevant?),
    # but the downside is it will make the compiled blast file huge (mitigate with eval).

    blastn \\
        -query ${query_fasta} \\
        -out ${query_name}.blast.out.raw \\
        -outfmt '6 ${BLAST_FORMAT}' \\
        -num_threads ${params.cpu_per_task} \\
        -db ref_db \\
        ${params.blast_args}
    
    # Filter and adding query_name to the first column
    awk -v qname="${query_name}" \\
        -v min_identity="${min_identity}" \\
        -v min_alignment="${min_alignment}" \\
        -v min_coverage="${min_coverage}" '
        BEGIN { OFS="\\t"; }
        {
            percent_identity = \$3;
            alignment_length = \$4;
            subject_length = \$17;
            subject_coverage = (alignment_length / subject_length) * 100;
            
            if (percent_identity >= min_identity && 
                alignment_length >= min_alignment && 
                subject_coverage >= min_coverage) {
                # Prepend query_name and print the entire line
                printf "%s\\t%s\\n", qname, \$0;
            }
        }' ${query_name}.blast.out.raw > ${query_name}.blast.out

    # Remove the raw output file if not needed
    if [ "${params.save_intermediate_files}" != "true" ]; then
        rm ${query_name}.blast.out.raw
    fi
    """
}

process RUN_KMA_SINGLE {
    publishDir "${params.results_directory}/2_kma", mode: 'copy'

    input:
    tuple val(sample_id), path(query_fastq), path("kma_db.comp.b"), path("kma_db.length.b"), path("kma_db.name"), path("kma_db.seq.b")

    output:
    path "${sample_id}.res", emit: kma_res
    path "${sample_id}.mat.gz", emit: kma_mat_gz

    script:
    """
    # Run KMA with single-end mode
    kma -i ${query_fastq} \\
        -o ${sample_id} \\
        -t_db kma_db \\
        -matrix \\
        -na -nc -nf \\
        -t ${params.cpu_per_task} \\
        ${params.kma_args}
    """
}

process RUN_KMA_PAIRED {
    publishDir "${params.results_directory}/2_kma", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq_files), path("kma_db.comp.b"), path("kma_db.length.b"), path("kma_db.name"), path("kma_db.seq.b")


    output:
    path "${sample_id}.res", emit: kma_res
    path "${sample_id}.mat.gz", emit: kma_mat_gz

    script:
    """
    # Run KMA with paired-end mode
    kma -ipe ${fastq_files[0]} ${fastq_files[1]} \\
        -o ${sample_id} \\
        -t_db kma_db \\
        -matrix \\
        -na -nc -nf \\
        -t ${params.cpu_per_task} \\
        ${params.kma_args}

    """
}

process COMPILE_BLAST_RESULTS {
    publishDir "${params.results_directory}/3_results", mode: 'copy'

    input:
    path blastout_files

    output:
    path "compiled_results.tsv", emit: compiled_results
    path "blast_coverage.tsv", emit: blast_coverage

    script:
    def headers = "query_file\t" + BLAST_FORMAT.trim().replaceAll("\\s+", "\t")
    """
    # Create unsorted compilation
    echo -e "${headers}" > compiled_results.tsv
    cat ${blastout_files} >> compiled_results.tsv

    # Process BLAST results to calculate coverage metrics
    processBLASTresults.py compiled_results.tsv blast_coverage.tsv
    """
}

process PRECLUSTER_QUERIES {
    publishDir "${params.results_directory}/2b_precluster_queries", mode: 'copy'

    input:
    path query_fastas
    path ref_fna

    output:
    path "query_features.tsv", emit: query_features
    path "clusters.tsv", emit: clusters
    path "representatives.tsv", emit: representatives
    path "cluster_manifest.tsv", emit: cluster_manifest
    path "non_aligning_queries.tsv", emit: non_aligning_queries
    path "hdbscan_tree.tsv", emit: hdbscan_tree

    script:
    def query_list = query_fastas.collect { it.getName() }.join('\n')
    def precluster_args = params.precluster_queries_args ? "--skani-args ${params.precluster_queries_args}" : ''
    """
    cat <<'EOF' > query_paths.txt
${query_list}
EOF

    preclusterQueries.py \\
        --reference ${ref_fna} \\
        --query-list query_paths.txt \\
        --threads ${params.cpu_per_task} \\
        ${precluster_args}
    """
}

process COMPILE_KMA_RESULTS {
    publishDir "${params.results_directory}/3_compiled_kma_results", mode: 'copy'

    input:
    path kma_res_files
    path kma_mat_files

    output:
    path "compiled_kma_res.tsv", emit: compiled_kma_res
    path "compiled_kma_mat.tsv", emit: compiled_kma_mat
    path "compiled_kma_mat_norm_1.tsv", emit: compiled_kma_mat_norm_1
    path "compiled_kma_mat_norm_2.tsv", emit: compiled_kma_mat_norm_2
    path "kma_coverage.tsv", emit: kma_coverage

    script:
    """
    # Create header
    head -n 1 \$(ls *.res | head -n 1) > header.tmp
    echo -e "SampleID\t\$(cat header.tmp)" > compiled_kma_res.tsv
    
    # Compile
    for res_file in *.res; do
        sample_id=\$(basename \$res_file .res)
        awk -v sid="\$sample_id" 'NR>1 {print sid "\\t" \$0}' \$res_file >> compiled_kma_res.tsv
    done
    
    # Process matrix files using the separate Python script
     processKmaMatrix.py --kma_outputs *.mat.gz \\
        --res_file compiled_kma_res.tsv \\
        --min_reads ${params.kma_min_reads} \\
        --max_reads ${params.kma_max_reads} \\
        --output compiled_kma_mat.tsv \\
        --output_normalised_1 compiled_kma_mat_norm_1.tsv \\
        --output_normalised_2 compiled_kma_mat_norm_2.tsv \\
        --coverage_output kma_coverage.tsv \\
        --window_size ${params.kma_window}
    """
}

process BUILD_CGVIEW_JSON_BLAST {
    publishDir "${params.results_directory}/3_results", mode: 'copy'

    input:
    path blast_file
    path ref_fasta
    path gff_file

    output:
    path "cgview.json"

    script:
    def gff_param = gff_file.name != 'NO_FILE' ? "--features ${gff_file}" : ''

    """
    cgviewBuilderPy.py \\
        --name "${params.cgview_title ?: 'CGViewMap'}" \\
        --sequence ${ref_fasta} \\
        --blast ${blast_file} \\
        ${gff_param} \\
        --output cgview.json
    """
}

process BUILD_CGVIEW_JSON_KMA {
    publishDir "${params.results_directory}/3_results", mode: 'copy'

    input:
    path kma_results_file
    path ref_fasta
    path gff_file

    output:
    path "cgview.json"

    script:
    def gff_param = gff_file.name != 'NO_FILE' ? "--features ${gff_file}" : ''

    """
    cgviewBuilderPy.py \\
        --name "${params.cgview_title ?: 'CGViewMap'}" \\
        --sequence ${ref_fasta} \\
        --plots ${kma_results_file} \\
        ${gff_param} \\
        --output cgview.json
    """
}

process GENERATE_SCIFR_REPORT {
    publishDir "${params.results_directory}/3_results", mode: "copy"

    input:
    path cgview_json
    path metadata
    path coverage_data
    path template

    output:
    path ("blitsfr.json"), optional: true, emit: scifr_input_json
    path ("blitsfr.html"), emit: scifr_report

    script:
    def metadata_param = metadata && metadata.name != '[]' ? "-m ${metadata}" : ""
    def coverage_param = coverage_data && coverage_data.name != '[]' ? "-d ${coverage_data}" : ""
    def save_intermediate_param = params.save_intermediate_files ? "--save_intermediate" : ""
    """
    echo '${groovy.json.JsonOutput.toJson(params)}' > params.json
    generateSCIFR.py \
        -c ${cgview_json} \
        ${metadata_param} \
        ${coverage_param} \
        ${save_intermediate_param} \
        -t ${template} \
        -o "blitsfr.html" \
        -j "blitsfr.json" \
        -l params.json \
        -p "${params.pipeline_version}"
    """
}

process VALIDATE_REPORT {
    publishDir "${params.results_directory}", mode: 'copy', pattern: "*.{html,html.gz}"
    
    input:
    path html_report
    
    output:
    path "blitsfr.result.html", emit: validated_html
    path "blitsfr.result.html.gz", emit: validated_html_gz
    
    script:
    """
    #!/usr/bin/env bash
    if awk '/@@BLITSFR@@INPUT@@START@@/ && /@@BLITSFR@@INPUT@@END@@/' ${html_report} | grep -q .; then
        echo "Validation passed: Report contains the expected JSON data block"
        cp -P ${html_report} blitsfr.result.html
        gzip -4 -c blitsfr.result.html > blitsfr.result.html.gz
        echo "Final report generated: blitsfr.result.html and blitsfr.result.html.gz"
    else
        echo "ERROR: VALIDATION FAILED - Report is not valid" >&2
        echo "The complete JSON data block with markers @@BLITSFR@@INPUT@@START@@ and @@BLITSFR@@INPUT@@END@@ was not found." >&2
        echo "This may be caused by error during template mutation, e.g. due to compute resource limitations." >&2
        echo "You may also want to verify there's sufficient memory available for the report generation step." >&2
        exit 1
    fi
    """
}

workflow.onComplete {
    def duration = workflow.duration
    def status = workflow.success ? 'SUCCESS' : 'FAILED'
    
    //ansi colour codes
    def green = '\033[32m'
    def red = '\033[31m'
    def blue = '\033[34m'
    def yellow = '\033[33m'
    def cyan = '\033[36m'
    def bold = '\033[1m'
    def reset = '\033[0m'
    
    def statusColour = workflow.success ? green : red
    
    println """
${blue}Started at:${reset} ${workflow.start}
${blue}Completed at:${reset} ${workflow.complete}
${yellow}Duration:${reset} ${bold}${duration}${reset}
${blue}Status:${reset} ${statusColour}${bold}${status}${reset}
"""
}
