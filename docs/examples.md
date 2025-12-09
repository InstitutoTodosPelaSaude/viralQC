# Practical Examples

## Example 1: Basic Dengue Analysis

```bash
# 1. Prepare environment (run once)
vqc get-nextclade-datasets
vqc get-blast-database

# 2. Run analysis
vqc run-from-fasta \
  --sequences-fasta dengue_samples.fasta \
  --output-dir qc_dengue \
  --output-file dengue_results.tsv \
  --cores 4
```

## Example 2: Influenza (Multiple Segments)

```bash
vqc run-from-fasta \
  --sequences-fasta influenza_samples.fasta \
  --output-dir qc_influenza \
  --output-file influenza_qc.tsv \
  --cores 8
```

ViralQC automatically identifies different segments (HA, NA, PB1, etc.).

## Example 3: Metagenomic Analysis

For unknown viruses, use relaxed parameters:

```bash
vqc run-from-fasta \
  --sequences-fasta metagenome.fasta \
  --output-dir qc_metagenome \
  --blast-pident 30 \
  --blast-qcov 30 \
  --blast-task dc-megablast \
  --ns-min-score 0.05
```

## Example 4: Reproducible Database

Create a BLAST database with a specific release date:

```bash
# Create database with sequences up to a specific date
vqc get-blast-database --release-date 2023-06-15

# Run analysis with this database
vqc run-from-fasta \
  --sequences-fasta samples.fasta \
  --blast-database datasets/blast.fasta
```

## Example 5: JSON Output

```bash
vqc run-from-fasta \
  --sequences-fasta samples.fasta \
  --output-file results.json \
  --cores 4
```

## Example 6: Sensitive BLAST Search

For divergent viruses:

```bash
vqc run-from-fasta \
  --sequences-fasta samples.fasta \
  --blast-task blastn \
  --blast-pident 70 \
  --blast-evalue 1e-5 \
  --cores 4
```
