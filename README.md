# BLITSFR
BLITSFR (BLAST Interactive Tracks in Single-File Report) is a Nextflow pipeline that compares the similarity of multiple sequences using BLAST or KMA and generates a single-file interactive report with circular genome visualisation.

## Example report
Try and download example of interactive BLITSFR report file from [https://scifr.fordelab.com/blitsfr](https://scifr.fordelab.com/blitsfr) or [https://nalarbp.github.io/blitsfr](https://nalarbp.github.io/blitsfr)


## Features
- **Interactive visualisation**: Circular genome plots with CGView
- **Single-file output**: Self-contained HTML report with embedded data
- **Metadata integration**: Optional sample metadata for enhanced visualisation
- **Scalable processing**: CLI, configurable CPU usage and parallel processing
- **Resume capability**: Nextflow workflow resumption for interrupted runs

## Installation

### Supported platforms
- macOS
- Linux
- Windows (via WSL2 only)

### Prerequisites
- Mamba (recommended) or Conda ([install here](https://conda-forge.org/download/))
- Git ([install here](https://git-scm.com/downloads))

`mamba` is preferred because it resolves environments faster than `conda`. `conda-lock` is recommended for reproducible installs, but not required.

### Windows (via WSL2)
- Use a Linux distribution in WSL2, such as Ubuntu.
- Run all installation and `blitsfr` commands inside the WSL2 terminal.

### Quick Install (Recommended)
This installer expects `git` and either `mamba` or `conda` to already be installed on your system.

```bash
curl -fsSL https://raw.githubusercontent.com/nalarbp/blitsfr/main/install.sh | bash
```

After installation:
```bash
conda activate blitsfr
blitsfr -h  # Confirm installation
```

### Manual Installation

<details>
<summary>Click to expand manual installation steps</summary>

1. Clone this repository:
```bash
git clone https://github.com/nalarbp/blitsfr.git
cd blitsfr
```

2. Create and activate the environment:
```bash
conda-lock install -n blitsfr conda-lock.yml
# Fallback: mamba env create -f environment.yml
# Fallback: conda env create -f environment.yml

conda activate blitsfr
```

3. Install BLITSFR:
```bash
pip install -e .
blitsfr -h  # Confirm installation
```

</details>

### How to update

1. Navigate to your BLITSFR directory and activate the environment:
```bash
cd blitsfr  # or wherever you installed BLITSFR
conda activate blitsfr
```

2. Pull the latest changes and update:
```bash
git pull origin main
pip install -e . --force-reinstall
```

## Basic usage
This repository contains example inputs in [sample/](sample/):
- `sample/assemblies_mode/` for assembly-mode examples
- `sample/reads_mode/` for read-mode examples

Try assemblies mode:

```bash
ls sample/assemblies_mode/ #to see example input files for assemblies mode
blitsfr assemblies -r sample/assemblies_mode/Reference.gbff -q 'sample/assemblies_mode/*.fna'
```

Try reads mode:

```bash
ls sample/reads_mode/ #to see example input files for reads mode
blitsfr reads --reads-mode paired -r sample/reads_mode/Reference.gbff -q 'sample/reads_mode/*_R{1,2}.fastq.gz'
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
blitsfr reads --reads-mode paired -r sample/reads_mode/Reference.gbff -q 'sample/reads_mode/*_R{1,2}.fastq.gz'
```

### Optional parameters
**Common options:**
- `-m, --metadata`: Metadata file (TSV) with sample information
- `-o, --output`: Output directory (default: results)
- `--title`: Title for CGView visualisation (default: 'CGView Map')
- `--cpu_per_task`: CPU cores per task (default: 2)
- `--resume`: Resume previous Nextflow run
- `--nf-args`: Additional Nextflow arguments passed to Nextflow, for example `--nf-args "-c custom.config -with-report report.html"`

**Assembly mode options:**
- `--blast-args`: Additional BLAST arguments (default: '-dust no -evalue 1E-20')
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
blitsfr assemblies -r sample/assemblies_mode/Reference.gbff -q 'sample/assemblies_mode/*.fna' \
  -m sample/assemblies_mode/metadata_vre_st78.txt \
  -o comparative_analysis
```

**Reads mode with metadata:**
```bash
blitsfr reads --reads-mode paired \
  -r sample/reads_mode/Reference.gbff \
  -q 'sample/reads_mode/*_R{1,2}.fastq.gz' \
  -m sample/reads_mode/metadata_vre_st78.txt \
  -o read_mapping_results
```

## Input file formats
**Reference file**: GenBank format (.gbk, .gbff) must containing the sequence of reference genome

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

The exact contents depend on whether you run `assemblies` or `reads`. The final validated reports are written to the output directory root. Nextflow working files are created in the directory where you launch `blitsfr`.

### Assemblies mode

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
├── blitsfr.result.html        #Final interactive report
└── blitsfr.result.html.gz     #Compressed final report
```

### Reads mode

```
results/
├── 1_reference/
│   ├── kma_db/
│   ├── ref_features.gff3
│   └── ref.fna
├── 2_kma/
│   ├── sample1.res
│   ├── sample1.mat.gz
│   ├── sample2.res
│   └── sample2.mat.gz
├── 3_compiled_kma_results/
│   ├── compiled_kma_res.tsv
│   ├── compiled_kma_mat.tsv
│   ├── compiled_kma_mat_norm_1.tsv
│   ├── compiled_kma_mat_norm_2.tsv
│   └── kma_coverage.tsv
├── 3_results/
│   ├── blitsfr.html
│   └── cgview.json
├── blitsfr.result.html
└── blitsfr.result.html.gz
```

### Nextflow runtime files

These are created in the directory where you run the command, not inside `results/`:

```
work/
.nextflow/
.nextflow.log
timeline.html
trace.txt
```

**Main output files**
- `blitsfr.result.html`: Main interactive report with embedded visualisations
- `blitsfr.result.html.gz`: Compressed version for sharing
- `3_results/cgview.json`: CGView configuration and data for visualisation
- `3_results/compiled_results.tsv`: Assembly-mode BLAST summary table
- `3_compiled_kma_results/compiled_kma_res.tsv`: Read-mode KMA summary table

## Citation
If you use BLITSFR in your research, please cite:

```
Interactive data analysis and reporting with SCIFR: A Single-File SPA Approach
Budi Permana, Thom P. Cuddihy, Brian M. Forde
bioRxiv 2025.07.10.664259; doi: https://doi.org/10.1101/2025.07.10.664259 
```

## License
This project is licensed under the Apache 2.0 - see the [LICENSE](LICENSE) for details.

## Documentation
Additional documentation and troubleshooting notes will be expanded in the repository wiki.

## Support
- **Issues**: Report bugs and request features on [GitHub Issues](https://github.com/nalarbp/blitsfr/issues)
- **Contact**: b.permana@uq.edu.au

## Acknowledgements
- CGViewer.js, please cite: Grant, J.R. and P. Stothard, CGView.js: a JavaScript package for visualizing small genomes. Journal of Open Source Software, 2026. 11(122): p. 9930.
- Nextflow devs and community
- ReactJS devs and community 
- Core JS libraries (Jotai.js, Nivo.js, AgGrid.js, D3.js) devs and community

---

**Version**: v0.1.2 
