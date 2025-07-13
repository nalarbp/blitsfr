# BLITSFR
BLITSFR (BLAST Interactive Tracks in Single-File Report) is a Nextflow pipeline that compares the similarity of multiple sequences using BLAST or KMA and generates a single-file interactive report with circular genome visualisation.

## Example report
Try and download example of interactive BLITSFR report file from [https://scifr.fordelab.com/blitsfr](https://scifr.fordelab.com/blitsfr)

## Features
- **Interactive visualisation**: Circular genome plots with CGView
- **Single-file output**: Self-contained HTML report with embedded data
- **Metadata integration**: Optional sample metadata for enhanced visualisation
- **Scalable processing**: CLI, configurable CPU usage and parallel processing
- **Resume capability**: Nextflow workflow resumption for interrupted runs

## Installation

### Prerequisites
- Conda or Mamba (https://conda-forge.org/download/)
- Git (https://git-scm.com/downloads)

### Quick install
1. Clone this repository:
```bash
git clone https://github.com/nalarbp/blitsfr.git
cd blitsfr
```

2. Create and activate the conda environment:
```bash
mamba env create -f environment.yml
mamba activate blitsfr
```

3. Install BLITSFR:
```bash
pip install -e .
blitsfr -h #to confirm its installed properly
```

### How to update

1. Navigate to your BLITSFR directory and activate the environment:
```bash
cd blitsfr
mamba activate blitsfr
```

2. Pull the latest changes and update from the repository:
```bash
git pull origin main
pip install -e . --force-reinstall
```

## Basic usage
This repository contains input file examples located in [sample/](sample/) directory for you to try blitsfr. Run the following command:

```bash
ls sample/ #to see example of required input files
blitsfr assemblies -r sample/Reference.gbff -q 'sample/*.fna'
```

### Required parameters
**For both modes:**
- `-r, --reference`: Reference sequence file in GenBank format
- `-q, --queries`: Query sequences (FASTA for assemblies, FASTQ for reads)

**Assembly mode specific:**
```bash
blitsfr assemblies -r reference.gbk -q 'assemblies/*.fasta'
```
**Read mode specific:**
```bash
blitsfr reads --reads-mode paired -r reference.gbk -q 'reads/*_R{1,2}.fastq.gz'
```

### Optional parameters
**Common options:**
- `-m, --metadata`: Metadata file (TSV) with sample information
- `-o, --output`: Output directory (default: results)
- `--title`: Title for CGView visualisation (default: 'CGView Map')
- `--cpu_per_task`: CPU cores per task (default: 2)
- `--resume`: Resume previous Nextflow run
- `-c, --config`: Nextflow configuration file
- `--nf-args`: Additional Nextflow arguments

**Assembly mode options:**
- `--blast-args`: Additional BLAST arguments (default: '-dust yes -evalue 1E-20')
- `--blast-filter-min-identity`: Minimum identity % for filtering (default: 80)
- `--blast-filter-min-coverage`: Minimum coverage % for filtering (default: 0)
- `--blast-filter-min-alignment`: Minimum alignment length in bp (default: 200)

**Read mode options:**
- `--reads-mode`: Single or paired-end reads (default: paired)
- `--window`: Window size for KMA score averaging (default: 1000)
- `--min-reads`: Minimum read count for normalisation (default: 0)
- `--max-reads`: Maximum read count for normalisation (default: 100)
- `--kma-args`: Additional KMA arguments

### Examples
**Assemblies comparison:**
```bash
blitsfr assemblies -r reference.gbk -q 'genomes/*.fasta' \
  --blast-filter-min-identity 90 \
  --blast-filter-min-coverage 50 \
  -o assembly_results
```

**With metadata:**
```bash
blitsfr assemblies -r reference.gbk -q 'samples/*.fna' \
  -m metadata_vre_st78.txt \
  -o comparative_analysis
```

## Input file formats
**Reference file**: GenBank format (.gbk, .gbff) containing the reference genome

**Query files**: 
- Assembly mode: FASTA format (.fasta, .fna, .fa)
- Read mode: FASTQ format (.fastq, .fq), compressed (.gz) supported

**Metadata file**: TSV format with mandatory 'id' column matching query file basenames (without file extension). Additional columns are user-specified and will be displayed in the interactive report.

Example metadata.tsv:
```
id	Geolocation2	Lineage_group	Assembly_acc	Strain
2_QLD_GCA_022046545.1	Australia_QLD	2	GCA_022046545.1	M95768
1_NSW_GCA_022046125.1	Australia_NSW	1	GCA_022046125.1	M87432
3_VIC_GCA_022046789.1	Australia_VIC	2	GCA_022046789.1	M98765
```

**Important**: The 'id' column must match the basename of your query files. For example, if your assembly file is `2_QLD_GCA_022046545.1.fna`, the corresponding id should be `2_QLD_GCA_022046545.1`.

## Output structure

```
results/
├── 1_reference/
│   ├── blast_db/
│   │   ├── ref_db.ndb
│   │   ├── ref_db.nhr
│   │   ├── ref_db.nin
│   │   ├── ref_db.njs
│   │   ├── ref_db.not
│   │   ├── ref_db.nsq
│   │   ├── ref_db.ntf
│   │   └── ref_db.nto
│   ├── ref_features.gff3
│   └── ref.fna
├── 2_blast/
│   ├── sample1.blast.out
│   └── sample2.blast.out
├── 3_results/
│   ├── blast_coverage.tsv
│   ├── blitsfr.html
│   ├── cgview.json
│   └── compiled_results.tsv
├── work/                      #Nextflow working directory
├── .nextflow/                 #Nextflow metadata
├── .nextflow.log              #Nextflow execution log
├── blitsfr.result.html        #Final interactive report
└── blitsfr.result.html.gz     #Compressed final report
```

**Main output files:**
- `blitsfr.result.html`: Main interactive report with embedded visualisations
- `blitsfr.result.html.gz`: Compressed version for sharing
- `3_results/compiled_results.tsv`: Tab-separated summary of all results
- `3_results/cgview.json`: CGView configuration and data for visualisation

## Citation
If you use BLITSFR in your research, please cite:

```
[Coming soon]
```

## License
This project is licensed under the Apache 2.0 - see the [LICENSE](LICENSE) for details.

## Documentation
[TODO] Add more detailed documentation on [docs/](docs/) dir.

## Support
- **Issues**: Report bugs and request features on [GitHub Issues](https://github.com/nalarbp/blitsfr/issues)
- **Contact**: b.permana@uq.edu.au

## Acknowledgements
- CGViewer.js
- Nextflow devs and community
- Kraken2, Bracken, and MetaPhlAn4 authors 
- NCBI and GTDB teams
- ReactJS devs and community 
- Core JS libraries (Jotai.js, Nivo.js, AgGrid.js, D3.js) devs and community

---

**Version**: v0.1.1 