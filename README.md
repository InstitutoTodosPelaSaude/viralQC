## Install

### Dependencies

```bash
micromamba env create -f env.yml
micromamba activate viralQC
```

### viralQC

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

### run-from-fasta

This command runs several steps to identify viruses represented in the input FASTA file and executes Nextclade for each identified virus/dataset. There are two modes for viral identification:

    - nextclade (default): Run `nextclade sort`
    - blast: Run `blastn`.

#### run-from-fasta (nextclade)

```bash
vqc run-from-fasta --sequences-fasta test_data/sequences.fasta
```

Some parameters can be specified:

- `--output-dir` — Output directory name. **Default:** `output`
- `--datasets-dir` — Path to the local Nextclade datasets directory. **Default:** `datasets`
- `--ns-min-score` — Minimum score used by the Nextclade `sort` command. **Default:** `0.1`
- `--ns-min-hits` — Minimum number of hits for Nextclade to consider a dataset. **Default:** `10`
- `--cores` — Number of threads used in `nextclade sort` and `nextclade run`. **Default:** `1`

The output directory has the following structure:

```
├── <datasets>                    # Output from nextclade sort; sequences for each dataset split into sequences.fa files.
├── datasets_selected.tsv         # Formatted nextclade sort output showing the mapping between input sequences and local datasets.
├── <virus/dataset>.nextclade.tsv # Nextclade run output for each identified virus, including clade assignments and QC metrics.
├── unmapped_sequences.txt        # Names of input sequences that were not mapped to any virus.
└── viruses.tsv                   # Nextclade sort output showing the mapping between input sequences and remote
```

#### run-from-fasta (nextclade)

For a run-from-fasta analysis using BLAST to identify the corresponding virus, a BLAST database must first be created:

```bash
vqc get-blast-database
```

Some parameters can be specified:

- `--output-dir` — Path to store the BLAST database.. **Default:** `datasets`
- `--datasets-dir` — Path of local nextclade datasets. **Default:** `datasets`
- `--cores` — Number of threads used in `nextclade sort` and `nextclade run`. **Default:** `1`


Then:

```bash
vqc run-from-fasta --sequences-fasta test_data/sequences.fasta --sort-mode blast
```

Some parameters can be specified:

- `--output-dir` — Output directory name. **Default:** `output`
- `--datasets-dir` — Path to the local BLAST database directory. **Default:** `datasets`
- `--cores` — Number of threads used in `nextclade sort` and `nextclade run`. **Default:** `1`

The output directory has the following structure:

```bash
├── datasets_selected.tsv            # Mapping between input sequences and local datasets.
├── <virus/dataset>.nextclade.tsv    # Nextclade run output for each identified virus, including clade assignments and QC metrics.
├── sequences.<virus/dataset>.fasta  # Input sequences split by virus/dataset.
├── unmapped_sequences.txt           # Names of input sequences that were not mapped to any virus.
└── viruses.tsv                      # BLAST output.
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
