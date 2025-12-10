# Commands and Usage

ViralQC provides three main commands through the command-line interface (`vqc`).

## get-nextclade-datasets

Downloads and configures Nextclade datasets locally.

```{important}
This command must be run **at least once** before using `run`.
```

### Usage

```bash
vqc get-nextclade-datasets --cores 2
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--datasets-dir` | String | `datasets` | Directory where datasets will be stored |
| `--cores` | Integer | `1` | Number of threads/cores to use |
| `--verbose` | Boolean | `False` | Show snakemake logs |

### Output Structure

```
datasets/
├── nextclade_data/
│   ├── denv1/
│   ├── denv2/
│   └── ...
├── external_datasets/
│   └── zikav/
└── external_datasets_minimizers.json
```

---

## get-blast-database

Creates a local BLAST database containing all viral genomes from NCBI RefSeq.

```{important}
This command must be run **at least once** before using `run`.
```

### Usage

```bash
vqc get-blast-database
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--output-dir` | String | `datasets` | Directory where the BLAST database will be stored |
| `--release-date` | String | `None` | Filter sequences by release date (YYYY-MM-DD). Only sequences released on or before this date will be included |
| `--cores` | Integer | `1` | Number of threads/cores to use |
| `--verbose` | Boolean | `False` | Show snakemake logs |

### Release Date Filtering

The `--release-date` parameter allows you to create a reproducible BLAST database by filtering sequences based on their NCBI release date:

```bash
# Create database with all sequences released up to June 15, 2023
vqc get-blast-database --release-date 2023-06-15
```

**Behavior:**
- When `--release-date` is provided:
  - Only sequences with `release_date <= specified_date` are included
  - The specified date is used as the database version identifier
- When not provided:
  - All available RefSeq sequences are included
  - Current date is used as the database version identifier

This is useful for:
- **Reproducibility**: Recreate the same database at different points in time
- **Auditing**: Track which sequences were available at a specific date
- **Comparative studies**: Analyze how results change with database updates

### Database Version

The database version is recorded in `blast.tsv` metadata file:
- Format: `ncbi-refseq-virus_YYYY-MM-DD`
- Uses the `--release-date` value if provided, otherwise the current date

### Output Structure

```
datasets/
├── blast.fasta          # Reference sequences
├── blast.fasta.ndb      # BLAST database files
├── blast.fasta.nhr
├── blast.fasta.nin
├── blast.fasta.nsq
├── blast.tsv            # Metadata with version info
└── blast_gff/           # GFF3 files for generic analysis
```

---

## run

Main analysis command. Identifies viruses, performs quality control, and extracts target regions.

### Usage

```bash
vqc run --input my_sequences.fasta
```

### Required Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `--input` | String | Path to the input FASTA file |

### Output Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--output-dir` | String | `output` | Directory for output files |
| `--output-file` | String | `results.tsv` | Results file (`.tsv`, `.csv`, or `.json`) |

### Dataset Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--datasets-dir` | String | `datasets` | Path to Nextclade datasets directory |

### Nextclade Sort Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--ns-min-score` | Float | `0.1` | Minimum score for valid match |
| `--ns-min-hits` | Integer | `10` | Minimum hits for valid match |

### BLAST Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--blast-database` | String | `datasets/blast.fasta` | Path to BLAST database |
| `--blast-database-metadata` | String | `datasets/blast.tsv` | Path to BLAST metadata |
| `--blast-pident` | Integer | `80` | Minimum percent identity (0-100) |
| `--blast-evalue` | Float | `1e-10` | Maximum E-value |
| `--blast-qcov` | Integer | `80` | Minimum query coverage (0-100) |
| `--blast-task` | String | `megablast` | BLAST task type |

### BLAST Task Types

The `--blast-task` parameter controls the BLAST algorithm sensitivity:

| Task | Description | Use Case |
|------|-------------|----------|
| `megablast` | Highly similar sequences (default) | Fast, same species |
| `dc-megablast` | Discontiguous megablast | Cross-species, more sensitive |
| `blastn` | Traditional BLASTN | More distant sequences |
| `blastn-short` | Short sequences | Sequences < 50 bp |

**Examples:**

```bash
# Default (megablast) - fast, for similar sequences
vqc run --input seqs.fasta

# More sensitive search for distant viruses
vqc run --input seqs.fasta --blast-task dc-megablast

# Traditional BLASTN for divergent sequences
vqc run --input seqs.fasta --blast-task blastn
```

### System Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--cores` | Integer | `1` | Number of threads/cores |
| `--verbose` | Boolean | `False` | Show snakemake logs |

### Complete Example

```bash
vqc run \
  --input samples.fasta \
  --output-dir results \
  --output-file report.tsv \
  --blast-pident 75 \
  --blast-task dc-megablast \
  --cores 8
```

### Analysis Workflow

1. **Nextclade Sort**: Maps sequences to local datasets
2. **BLAST Analysis**: Identifies unmapped sequences
3. **Nextclade Run**: Quality control analysis
4. **Post-processing**: Combines and scores results
5. **Region Extraction**: Extracts target regions based on quality

## API

How the viralQC is designed to integrate with viral genomic databases, it is possible to integrate the analysis module in the code of other applications.

This can be done by importing the `RunAnalysis` class from the `viralqc.core.run_analysis` module. This class has the `run` method that executes the quality analysis of a viral genome, receiving as parameter the path to the FASTA file containing the sequences to be analyzed. Other parameters can be informed in an optimized way.

### Usage

!important  

```python
from viralqc.core.run_analysis import RunAnalysis

input_file = "seqs.fasta"
output_directory = "results"
output_file = "results.json"
run_analysis = RunAnalysis()

snakemake_response = run_analysis.run(
        sequences_fasta=input_file,
        output_dir=output_directory,
        output_file=output_file
)
```

Or a flexible approach:

```python
from viralqc.core.run_analysis import RunAnalysis

input_file = "seqs.fasta"
output_directory = "results"
output_file = "results.json"
run_analysis = RunAnalysis()

snakemake_response = run_analysis.run(
        sequences_fasta=input_file,
        output_dir=output_directory,
        output_file=output_file,
        cores=2,
        datasets_local_path="datasets",
        nextclade_sort_min_score=0.1,
        nextclade_sort_min_hits=10,
        blast_database="datasets/blast.fasta",
        blast_database_metadata="datasets/blast.tsv",
        blast_identity_threshold=0,
        blast_evalue=0.01,
        blast_qcov=0,
        blast_task="blastn"
)
```

To check the results:

```python
if snakemake_response.status == 200:
    results_data = snakemake_response.get_results()
    for seq_result in results_data:
        virus = seq_result.get("virus")
        quality = seq_result.get("genomeQuality")
        coverage = seq_result.get("coverage")
        print(virus, quality, coverage)
else:
    raise Exception(snakemake_response.format_log())
```

### Attributes and Methods

The `run` method returns a `SnakemakeResponse` object that has the following attributes:

| Attribute | Type | Description |
|-----------|------|-------------|
| run_id | str | Execution ID |
| status | RunStatus | Execution status, which can be 200 (success) or 500 (failure) |
| log_path | str | Path to the log file |
| results_path | str | Path to the results file |

And the following methods:

| Method | Description |
|--------|-----------|
| `format_log()` | Returns the log file content formatted |
| `get_results()` | Returns the results file content in dictionary format |
