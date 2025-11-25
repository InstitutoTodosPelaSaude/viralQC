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
  "seqtk>=1.5.0,<1.6.0" \
  "ncbi-datasets-cli>=18.9.0,<18.10.0" \
  "taxonkit>=0.20.0,<0.21.0"
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
vqc get-blast-database --cores 2 --output-dir <directory_name>
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
output/
├── identified_datasets/              # Outputs from nextclade sort
│   ├── datasets_selected.tsv         # Formatted nextclade sort output showing the mapping between input sequences and local datasets
│   ├── viruses.tsv                   # Nextclade sort output showing the mapping between input sequences and remote (nextclade_data) datasets
│   ├── viruses.external_datasets.tsv # Nextclade sort output showing the mapping between input sequences and external (outside nextclade_data) datasets
│   ├── unmapped_sequences.txt        # Names of input sequences that were not mapped to any virus on nextclade sort
│   └── <virus/dataset>/              # Sequences for each dataset split into sequences.fa files
│       └── sequences.fa
├── blast_results/                    # BLAST analysis results
│   ├── unmapped_sequences.blast.tsv  # BLAST results for unmapped sequences
│   └── blast_viruses.list            # List of unique virus accessions identified by BLAST
├── nextclade_results/                # Nextclade run outputs
│   ├── <virus>.nextclade.tsv         # Standard nextclade run output for each identified virus
│   └── <accession>.generic.nextclade.tsv # Generic nextclade run output for BLAST-identified viruses
├── gff_files/                        # GFF annotation files
│   ├── <virus>.nextclade.gff         # GFF files from standard nextclade runs
│   └── <accession>.generic.nextclade.gff # GFF files from generic nextclade runs
├── logs/                             # Log files for all pipeline steps
│   ├── nextclade_sort.log
│   ├── blast.log
│   ├── nextclade.<virus>.log
│   ├── generic_nextclade.<virus>.log
│   ├── pp_nextclade.log
│   └── extract_target_regions.log
├── results(.tsv,.csv,.json)                       # Final combined results (format: .tsv, .csv, or .json)
├── sequences_target_regions.bed      # BED file with target regions (includes 4 columns: seqname, start, end, region_name)
└── sequences_target_regions.fasta    # FASTA file with extracted target region sequences
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

## Commands Logic

### get-nextclade-datasets

This command is responsible for downloading and configuring the Nextclade datasets required for analysis. It reads a configuration file (`config/datasets.yml`) to determine which datasets to fetch.

![get-nextclade-datasets](assets/get_nextclade_datasets.svg)

**Logic:**
1.  **Read Configuration:** The command reads the `config/datasets.yml` file to identify the list of viruses and their corresponding Nextclade dataset names or GitHub repositories.
2.  **Download Nextclade Datasets:** For viruses hosted on Nextclade's official repository, it uses the `nextclade dataset get` command to download the reference sequence, gene map, and other necessary files.
3.  **Download GitHub Datasets:** For viruses hosted on GitHub, it uses a custom Python script (`get_github_dataset.py`) to download the dataset files.
4.  **Generate Minimizer Index:** For the GitHub datasets, it generates a minimizer index (`external_datasets_minimizers.json`) using `get_minimizer_index.py`. This index is crucial for the `nextclade sort` step to recognize these external datasets.

### get-blast-database

This command builds a local BLAST database containing all viral genomes from NCBI RefSeq. This database is used to identify sequences that are not recognized by Nextclade's default datasets.

![get-blast-database](assets/get_blast_database.svg)

**Logic:**
1.  **Download Genomes:** It downloads all viral genomes from NCBI RefSeq using the `datasets` CLI tool.
2.  **Format Database:** It processes the downloaded files to create a FASTA file (`blast.fasta`) and extracts metadata (accession, virus name, taxonomy ID) into a TSV file (`blast.tsv`).
3.  **Taxonomy Information:** It downloads NCBI taxonomy data and uses `taxonkit` to add species names and taxonomy IDs to the metadata.
4.  **Create BLAST DB:** It uses `makeblastdb` to create the searchable BLAST database.
5.  **Generate GFFs:** It converts the JSONL annotation reports from NCBI into GFF3 format using `jsonl_to_gff.py`. **Crucially**, this step filters out accessions where the CDS length is not divisible by 3 or other structural anomalies, ensuring that only valid annotations are used for generic Nextclade runs.

### run-from-fasta

This is the main analysis command. It takes a FASTA file of query sequences and performs virus identification, clade assignment, and quality control.

![run-from-fasta](assets/run_from_fasta.svg)

**Logic:**
1.  **Nextclade Sort:** The input sequences are first passed to `nextclade sort`. This tool compares the sequences against the local Nextclade datasets (including the external ones via the minimizer index) to identify which virus each sequence belongs to.
2.  **Datasets Identification:**: Identify the virus/dataset that each sequence belongs to with format_nextclade_sort.py
    *   **Mapped Sequences:** Sequences that match a known dataset are assigned to that virus.
    *   **Unmapped Sequences:** Sequences that do not match any dataset are flagged as "unmapped".
3.  **BLAST Analysis:** The unmapped sequences are queried against the local BLAST database (`blastn`) to identify their closest viral match.
4.  **Nextclade Analysis:**
    *   **A. Standard Run:** For sequences mapped to a known Nextclade dataset, `nextclade run` is executed using that dataset.
    *   **B. Generic Run:** For sequences identified by BLAST (but not in the Nextclade datasets), a "generic" `nextclade run` is performed. This uses the reference sequence and GFF annotation from the BLAST database (created by `get-blast-database`).
5.  **Post-Processing:** The results from all `nextclade run` executions (standard and generic) and the BLAST results for any remaining unmapped sequences are combined into a single output file (`results.tsv`).
6.  **Target Region Extraction:** Finally, specific genomic regions of interest (if defined) are extracted from the sequences based on quality criteria.

