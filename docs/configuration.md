# Dataset Configuration

ViralQC uses a configuration file called `datasets.yml` to define which viruses and datasets are available for analysis. This file is located at `viralqc/config/datasets.yml`.

## datasets.yml File Structure

The file has two main sections:

1. **`nextclade_data`**: Datasets hosted in the official Nextclade repository
2. **`github`**: Custom datasets hosted on GitHub

```yaml
nextclade_data:
  virus-identifier:
    dataset: "dataset/path"
    tag: "dataset-version"
    virus_name: "Virus Name"
    # ... other parameters

github:
  virus-identifier:
    repository: "user/repository"
    dataset: "dataset/path"
    tag: "development"
    virus_name: "Virus Name"
    # ... other parameters
```

## Nextclade Datasets

Nextclade datasets are official or community datasets available through the [nextclade_data](https://github.com/nextstrain/nextclade_data) repository.

```yaml
nextclade_data:
  denv1:
    dataset: "community/v-gen-lab/dengue/denv1"
    tag: "2025-04-02--19-11-08Z"
    virus_name: "Dengue virus type 1"
    virus_tax_id: 11053
    virus_species: "Orthoflavivirus denguei"
    virus_species_tax_id: 3052464
    segment: "Unsegmented"
    ncbi_id: "NC_001477.1"
    target_gene: "E"
    target_regions: ["C", "prM", "E"]
    private_mutation_total_threshold: 70
```

## GitHub Datasets

Custom datasets can be hosted in GitHub repositories.

```yaml
github:
  zikav:
    repository: "dezordi/nextclade_data_workflows"
    dataset: "zikaV/dataset"
    tag: "development"
    virus_name: "Zika virus"
    virus_tax_id: 64320
    virus_species: "Orthoflavivirus zikaense"
    virus_species_tax_id: 3048459
    segment: "Unsegmented"
    ncbi_id: "NC_035889.1"
    target_gene: "E"
    target_regions: ["C", "prM", "E"]
    private_mutation_total_threshold: 40
```

## Configuration Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `dataset` | String | Dataset path in nextclade_data or GitHub repositories |
| `repository` | String | GitHub repository name in `user/repo` format (for `github`) |
| `tag` | String | Dataset version/tag or repository branch |
| `virus_name` | String | Full virus name |
| `virus_tax_id` | Integer | Virus taxonomic ID in NCBI Taxonomy |
| `virus_species` | String | Viral species name |
| `virus_species_tax_id` | Integer | Species taxonomic ID |
| `segment` | String | Segment name (use "Unsegmented" for non-segmented viruses) |
| `ncbi_id` | String | Reference genome accession in NCBI |
| `target_gene` | String | Target gene/CDS name |
| `target_regions` | List | List of target genes/CDS |
| `private_mutation_total_threshold` | Integer | Private mutation threshold for quality control |

```{note}
For non-segmented viruses (e.g., Dengue, Zika, SARS-CoV-2), use `"Unsegmented"` for the segment field.
For segmented viruses, specify the segment name (e.g., `"HA"`, `"NA"`, `"L"`, `"M"`, `"S"`).
```
