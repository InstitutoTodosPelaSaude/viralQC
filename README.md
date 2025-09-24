## Install

ViralQC is a tool and package for quality control of consensus virus genomes which uses the [nextclade tool](https://docs.nextstrain.org/projects/nextclade/en/stable/), [BLAST](https://www.ncbi.nlm.nih.gov/books/NBK279690/) and a series of internal logics to classify viral sequences and perform quality control of complete genomes, regions or target genes.

### From pip

First, install the dependencies:

- [Snakemake 7.32](https://snakemake.readthedocs.io/en/v7.32.0/getting_started/installation.html)
- [Ncbi-blast 2.16](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.16.0/)
- [Nextclade 3.15](https://docs.nextstrain.org/projects/nextclade/en/3.15.3/user/nextclade-cli/installation/)
- [Seqtk 1.5](https://github.com/lh3/seqtk/releases/tag/v1.5)
- [Python < 3.12](https://www.python.org/downloads/)

or with micromamba

```bash
micromamba install \
  -c conda-forge \
  -c bioconda \
  "python>=3.8.0,<3.12.0" \
  "snakemake-minimal>=7.32.0,<7.33.0" \
  "blast>=2.16.0,<2.17.0" \
  "nextclade>=3.15.0,<3.16.0" \
  "seqtk>=1.5.0,<1.6.0"
```

Then, install viralQC

```bash
pip install viralQC
```

### From Source

```bash
git clone https://github.com/InstitutoTodosPelaSaude/viralQC.git
cd viralQC
```

#### Dependencies

```bash
micromamba env create -f env.yml
micromamba activate viralQC
```

#### viralQC

```bash
pip install .
```

### Check installation (CLI)

```bash
vqc --help
```

## Usage (CLI)

### get-nextclade-datasets

This command configures local datasets using nextclade. It is necessary to run at least once to generate a local copy of the nextclade datasets, before running the `run-from-fasta` command

```bash
vqc get-nextclade-datasets --cores 2
```

A directory name can be specified, the default is `datasets`.

```bash
vqc get-nextclade-datasets --cores 2 --datasets-dir <directory_name>
```

### get-blast-database

This command configures local blast database with all ncbi refseq viral genomes. It is necessary to run at least once to generate a local blast database, before running the `run-from-fasta` command.

```bash
vqc get-blast-database --cores 2
```

A output directory name can be specified, the default is `datasets`.

```bash
vqc get-nextclade-datasets --cores 2 --output-dir <directory_name>
```

### run-from-fasta

This command runs several steps to identify viruses represented in the input FASTA file and executes Nextclade for each identified virus/dataset.

#### run-from-fasta

```bash
vqc run-from-fasta --sequences-fasta test_data/sequences.fasta
```

Some parameters can be specified:

- `--output-dir` — Output directory name. **Default:** `output`
- `--output-file` - File to write final results. Valid extensions: .csv, .tsv or .json. **Default:** `results.tsv`
- `--datasets-dir` — Path to the local Nextclade datasets directory. **Default:** `datasets`
- `--ns-min-score` — Minimum score used by the Nextclade `sort` command. **Default:** `0.1`
- `--ns-min-hits` — Minimum number of hits for Nextclade to consider a dataset. **Default:** `10`
- `--blast-database` - Path to store local blast database. **Default:** `datasets/blast.fasta`
- `--identity-threshold` - Percentual identity threshold for BLAST analysis. **Default:** `0.9`
- `--cores` — Number of threads used in `nextclade sort` and `nextclade run`. **Default:** `1`

The output directory has the following structure:

```
├── <datasets>                    # Output from nextclade sort; sequences for each dataset split into sequences.fa files.
├── datasets_selected.tsv         # Formatted nextclade sort output showing the mapping between input sequences and local datasets.
├── <virus/dataset>.nextclade.tsv # Nextclade run output for each identified virus, including clade assignments and QC metrics.
├── unmapped_sequences.txt        # Names of input sequences that were not mapped to any virus on nextclade sort.
├── unmapped_sequences.blast.tsv  # BLAST results for unmapped sequences.
└── viruses.tsv                   # Nextclade sort output showing the mapping between input sequences and remote
```

## Usage (API)

```bash
vqc-server
```

Go to `http://127.0.0.1:8000/docs`

### Development

Install development dependencies and run `black` into `viralqc` directory.

```bash
pip install -e ".[dev]"
black viralqc
```
