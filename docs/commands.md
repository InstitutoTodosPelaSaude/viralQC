# Commands and Usage

ViralQC provides three main commands through the command-line interface (`vqc`).

## get-nextclade-datasets

Downloads and configures Nextclade datasets locally.

```{important}
This command must be run **at least once** before using `run-from-fasta`.
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
This command must be run **at least once** before using `run-from-fasta`.
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

## run-from-fasta

Main analysis command. Identifies viruses, performs quality control, and extracts target regions.

### Usage

```bash
vqc run-from-fasta --sequences-fasta my_sequences.fasta
```

### Required Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `--sequences-fasta` | String | Path to the input FASTA file |

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
vqc run-from-fasta --sequences-fasta seqs.fasta

# More sensitive search for distant viruses
vqc run-from-fasta --sequences-fasta seqs.fasta --blast-task dc-megablast

# Traditional BLASTN for divergent sequences
vqc run-from-fasta --sequences-fasta seqs.fasta --blast-task blastn
```

### System Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--cores` | Integer | `1` | Number of threads/cores |

### Complete Example

```bash
vqc run-from-fasta \
  --sequences-fasta samples.fasta \
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
