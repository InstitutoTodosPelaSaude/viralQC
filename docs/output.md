# Output Structure

When you run `run-from-fasta`, ViralQC creates the following output:

```
output/
├── identified_datasets/
│   ├── datasets_selected.tsv
│   ├── viruses.tsv
│   ├── viruses.external_datasets.tsv
│   ├── unmapped_sequences.txt
│   └── <virus>/sequences.fa
├── blast_results/
│   ├── unmapped_sequences.blast.tsv
│   └── blast_viruses.list
├── nextclade_results/
│   ├── <virus>.nextclade.tsv
│   └── <accession>.generic.nextclade.tsv
├── gff_files/
│   ├── <virus>.nextclade.gff
│   └── <accession>.generic.nextclade.gff
├── logs/
│   ├── nextclade_sort.log
│   ├── blast.log
│   └── ...
├── results.tsv
├── sequences_target_regions.bed
└── sequences_target_regions.fasta
```

## Main Results File

The `results.tsv` contains:

### Identification Columns

| Column | Description |
|--------|-------------|
| `seqName` | Sequence name from FASTA |
| `virus` | Identified virus name |
| `virus_tax_id` | NCBI taxonomic ID |
| `segment` | Genomic segment |
| `clade` | Phylogenetic clade |
| `datasetVersion` | Dataset version used |

### Quality Columns (ViralQC)

| Column | Description |
|--------|-------------|
| `genomeQuality` | Overall quality (A, B, C, D) |
| `genomeQualityScore` | Numeric score (0-24) |
| `missingDataQuality` | Missing data quality |
| `privateMutationsQuality` | Private mutations quality |
| `cdsCoverageQuality` | CDS coverage quality |
| `targetRegionsQuality` | Target regions quality |

### Mutation Columns

| Column | Description |
|--------|-------------|
| `totalSubstitutions` | Total nucleotide substitutions |
| `totalDeletions` | Total deletions |
| `totalInsertions` | Total insertions |
| `coverage` | Genome coverage (0.0-1.0) |

## Target Regions Files

### sequences_target_regions.bed

```
seq1    94      2419    C,prM,E
seq2    0       10735   genome
```

### sequences_target_regions.fasta

Extracted sequences from regions meeting quality criteria.
