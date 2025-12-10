# How to Add New Datasets

## Adding a Nextclade Dataset

If the dataset is available in the official Nextclade repository:

### Step 1: Identify the Dataset

```bash
nextclade dataset list
```

### Step 2: Edit datasets.yml

Add a new entry in `viralqc/config/datasets.yml`:

```yaml
nextclade_data:
  my-new-virus:
    dataset: "complete/dataset/path"
    tag: "2025-XX-XX--XX-XX-XXZ"
    virus_name: "Full Virus Name"
    virus_tax_id: 123456
    virus_species: "Species Name"
    virus_species_tax_id: 789012
    segment: "Unsegmented"
    ncbi_id: "NC_XXXXXX.X"
    target_gene: "gene_name"
    target_regions: ["gene1", "gene2"]
    private_mutation_total_threshold: 50
```

### Step 3: Obtain Taxonomic Information

Consult [NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy) for:
- `virus_tax_id`: Virus taxonomic ID
- `virus_species_tax_id`: Species taxonomic ID

### Step 4: Test

```bash
vqc get-nextclade-datasets --cores 2
```

---

## Adding a GitHub Dataset

For custom datasets hosted on GitHub:

### Step 1: Prepare Repository

Your repository must contain:

```
your-repository/
└── dataset/path/
    ├── reference.fasta
    ├── genome_annotation.gff3
    ├── tree.json (optional)
    ├── pathogen.json
    └── sequences.fasta (optional)
```

### Step 2: Add to datasets.yml

```yaml
github:
  my-custom-virus:
    repository: "your-user/your-repository"
    dataset: "path/within/repo"
    tag: "main"
    virus_name: "My Custom Virus"
    virus_tax_id: 123456
    virus_species: "Species Name"
    virus_species_tax_id: 789012
    segment: "Unsegmented"
    ncbi_id: "NC_XXXXXX.X"
    target_gene: "gene1"
    target_regions: ["gene1", "gene2"]
    private_mutation_total_threshold: 40
```

### Step 3: Test

```bash
vqc get-nextclade-datasets --cores 2
```

Verify dataset downloaded to `datasets/external_datasets/my-custom-virus/`.
