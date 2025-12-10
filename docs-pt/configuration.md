# Configuração de Datasets

O ViralQC utiliza um arquivo de configuração chamado `datasets.yml` para definir quais vírus e datasets estão disponíveis. Este arquivo está localizado em `viralqc/config/datasets.yml`.

## Estrutura do arquivo datasets.yml

O arquivo possui duas seções principais:

1. **`nextclade_data`**: Datasets do repositório oficial do Nextclade
2. **`github`**: Datasets personalizados hospedados no GitHub

```yaml
nextclade_data:
  virus-identificador:
    dataset: "caminho/do/dataset"
    tag: "versão-do-dataset"
    virus_name: "Nome do Vírus"
    # ... outros parâmetros

github:
  virus-identificador:
    repository: "usuario/repositorio"
    dataset: "caminho/do/dataset"
    tag: "development"
    virus_name: "Nome do Vírus"
    # ... outros parâmetros
```

## Datasets do Nextclade

Datasets oficiais ou da comunidade disponíveis no [nextclade_data](https://github.com/nextstrain/nextclade_data):

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

## Datasets do GitHub

Datasets personalizados hospedados em repositórios GitHub:

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

## Parâmetros de Configuração

| Parâmetro | Tipo | Descrição |
|-----------|------|-----------|
| `dataset` | String | Caminho do dataset no nextclade_data ou GitHub |
| `repository` | String | Nome do repositório GitHub (`usuario/repo`) |
| `tag` | String | Versão/tag do dataset ou branch |
| `virus_name` | String | Nome completo do vírus |
| `virus_tax_id` | Integer | ID taxonômico do vírus no NCBI |
| `virus_species` | String | Nome da espécie viral |
| `virus_species_tax_id` | Integer | ID taxonômico da espécie |
| `segment` | String | Nome do segmento ("Unsegmented" para não-segmentados) |
| `ncbi_id` | String | Accession do genoma de referência |
| `target_gene` | String | Nome do gene/CDS alvo |
| `target_regions` | List | Lista de genes/CDS alvo |
| `private_mutation_total_threshold` | Integer | Limite de mutações privadas |

```{note}
Para vírus não-segmentados (Dengue, Zika, SARS-CoV-2), use `"Unsegmented"`.
Para vírus segmentados, especifique o nome do segmento (ex: `"HA"`, `"NA"`, `"L"`, `"M"`, `"S"`).
```
