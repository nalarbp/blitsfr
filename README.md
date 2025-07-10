# BLITSFR
BLITSFR (BLAST Interactive Tracks in Single-File Report) is a Nextflow pipeline that compare the similarity of multiple sequences using BLAST and generate a single-file interactive report.

## Example report
Try and download example of interactive BLITSFR report file from [https://scifr.fordelab.com/blitsfr](https://scifr.fordelab.com/blitsfr)

## Features

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

## Basic usage
This repository contains input file examples located in [sample/](sample/) directory for you to try blitsfr. Run the following command:

```bash
ls sample/ #to see example of required input files
blitsfr assemblies -r sample/Reference.gbff -q 'sample/*.fna'
```

### Required parameters


### Optional parameters


### Examples

## Input file formats

Example of input files are available in [sample/](sample/) directory.

## Output structure

```
results/
...
├── blitsfr.result.html #Final report
└── blitsfr.result.html.gz #Compressed final report
```

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

**Version**: v0.1.0 