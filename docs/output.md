# Output Structure

When you run `run-from-fasta`, ViralQC creates the following output:

```
outputs/             # User specified output directory (e.g., --output-dir my_results)
└── .snakemake/     # Snakemake run files
└── outputs/        # ViralQC output files
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

### Main File: results.tsv (or .csv, .json)

This is the file containing consolidated results from all analyses:

##### 1. Sequence Identification


| Column | Type | Description |
|--------|------|-------------|
| `seqName` | String | Sequence name in the input FASTA file |
| `virus` | String | Identified virus name |
| `virus_tax_id` | Integer | Virus taxonomic ID in NCBI Taxonomy |
| `virus_species` | String | Viral species name |
| `virus_species_tax_id` | Integer | Species taxonomic ID |
| `segment` | String | Genomic segment (e.g., "HA", "NA", "Unsegmented") |
| `ncbi_id` | String | Reference genome accession in NCBI |
| `dataset` | String | Dataset identifier used |
| `datasetVersion` | String | Dataset version/tag |
| `clade` | String | Phylogenetic clade (when available) |


##### 2. Quality Metrics (ViralQC A-D System)

| Column | Type | Description |
|--------|------|-------------|
| `genomeQuality` | String | Overall genome quality (A, B, C, or D) |
| `genomeQualityScore` | Float | Normalized numeric score (0-24) |
| `missingDataQuality` | String | Quality based on missing data (A, B, C, or D) |
| `privateMutationsQuality` | String | Quality based on private mutations (A, B, C, or D) |
| `mixedSitesQuality` | String | Quality based on mixed sites (A, B, C, or D) |
| `snpClustersQuality` | String | Quality based on SNP clusters (A, B, C, or D) |
| `frameShiftsQuality` | String | Quality based on frameshifts (A, B, C, or D) |
| `stopCodonsQuality` | String | Quality based on stop codons (A, B, C, or D) |
| `cdsCoverageQuality` | String | Coverage quality per CDS (e.g., "E: A, prM: B") |
| `targetRegionsQuality` | String | Target regions quality (A, B, C, D, or empty) |
| `targetGeneQuality` | String | Target gene quality (A, B, C, D, or empty) |


##### 3. Quality Metrics (Nextclade)

| Column | Type | Description |
|--------|------|-------------|
| `qc.overallScore` | Float | Nextclade overall quality score |
| `qc.overallStatus` | String | Nextclade quality status (good, mediocre, bad) |
| `qc.privateMutations.total` | Integer | Total private mutations (Nextclade) |
| `qc.privateMutations.score` | Float | Private mutations score (Nextclade) |
| `qc.privateMutations.status` | String | Private mutations status (Nextclade) |
| `qc.missingData.score` | Float | Missing data score (Nextclade) |
| `qc.missingData.status` | String | Missing data status (Nextclade) |
| `qc.mixedSites.totalMixedSites` | Integer | Total mixed sites (Nextclade) |
| `qc.mixedSites.score` | Float | Mixed sites score (Nextclade) |
| `qc.mixedSites.status` | String | Mixed sites status (Nextclade) |
| `qc.snpClusters.totalSNPs` | Integer | Total clustered SNPs (Nextclade) |
| `qc.snpClusters.score` | Float | SNP clusters score (Nextclade) |
| `qc.snpClusters.status` | String | SNP clusters status (Nextclade) |
| `qc.frameShifts.totalFrameShifts` | Integer | Total frameshifts (Nextclade) |
| `qc.frameShifts.score` | Float | Frameshifts score (Nextclade) |
| `qc.frameShifts.status` | String | Frameshifts status (Nextclade) |
| `qc.stopCodons.totalStopCodons` | Integer | Total stop codons (Nextclade) |
| `qc.stopCodons.score` | Float | Stop codons score (Nextclade) |
| `qc.stopCodons.status` | String | Stop codons status (Nextclade) |

##### 4. Coverage and Regions

| Column | Type | Description |
|--------|------|-------------|
| `coverage` | Float | Genome coverage (0.0 to 1.0) |
| `cdsCoverage` | String | Coverage of each CDS (format: "gene1: 0.98, gene2: 1.0") |
| `targetRegionsCoverage` | String | Coverage of target regions defined in `target_regions` |
| `targetGeneCoverage` | String | Coverage of target gene defined in `target_gene` |
| `targetRegions` | String | List of target regions (separated by \|) |
| `targetGene` | String | Main target gene name |

##### 5. Nucleotide Mutations

| Column | Type | Description |
|--------|------|-------------|
| `totalSubstitutions` | Integer | Total nucleotide substitutions |
| `totalDeletions` | Integer | Total nucleotide deletions |
| `totalInsertions` | Integer | Total nucleotide insertions |
| `totalFrameShifts` | Integer | Total frameshift mutations |
| `totalMissing` | Integer | Total missing nucleotides (N's or gaps) |
| `totalNonACGTNs` | Integer | Total non-ACGTN characters |
| `substitutions` | String | List of substitutions (format: gene:pos:ref>alt) |
| `deletions` | String | List of deletions |
| `insertions` | String | List of insertions |
| `frameShifts` | String | List of frameshifts |
| `alignmentScore` | Float | Alignment score |

##### 6. Amino Acid Mutations

| Column | Type | Description |
|--------|------|-------------|
| `totalAminoacidSubstitutions` | Integer | Total amino acid substitutions |
| `totalAminoacidDeletions` | Integer | Total amino acid deletions |
| `totalAminoacidInsertions` | Integer | Total amino acid insertions |
| `totalUnknownAa` | Integer | Total unknown amino acids |
| `aaSubstitutions` | String | List of amino acid substitutions |
| `aaDeletions` | String | List of amino acid deletions |
| `aaInsertions` | String | List of amino acid insertions |

##### 7. Private Mutations (Detailed)

| Column | Type | Description |
|--------|------|-------------|
| `privateNucMutations.totalPrivateSubstitutions` | Integer | Total private substitutions |
| `privateNucMutations.totalLabeledSubstitutions` | Integer | Total known/cataloged private mutations |
| `privateNucMutations.totalUnlabeledSubstitutions` | Integer | Total uncataloged private mutations |
| `privateNucMutations.totalReversionSubstitutions` | Integer | Total reversions (mutations that revert to ancestral reference) |

**Note on output formats:**
- **TSV/CSV**: All columns are strings or numeric values
- **JSON**: Columns like `cdsCoverage`, `cdsCoverageQuality`, and `targetRegionsCoverage` are formatted as arrays of objects for easier programmatic parsing

## Target Regions Files

### sequences_target_regions.bed

```
seq1    94      2419    C,prM,E
seq2    0       10735   genome
```

### sequences_target_regions.fasta

Extracted sequences from regions meeting quality criteria.
