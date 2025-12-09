# Como Adicionar Novos Datasets

## Adicionando Dataset do Nextclade

Se o dataset está disponível no repositório oficial do Nextclade:

### Passo 1: Identificar o Dataset

```bash
nextclade dataset list
```

### Passo 2: Editar datasets.yml

Adicione uma nova entrada em `viralqc/config/datasets.yml`:

```yaml
nextclade_data:
  meu-novo-virus:
    dataset: "caminho/completo/do/dataset"
    tag: "2025-XX-XX--XX-XX-XXZ"
    virus_name: "Nome Completo do Vírus"
    virus_tax_id: 123456
    virus_species: "Nome da Espécie"
    virus_species_tax_id: 789012
    segment: "Unsegmented"
    ncbi_id: "NC_XXXXXX.X"
    target_gene: "nome_do_gene"
    target_regions: ["gene1", "gene2"]
    private_mutation_total_threshold: 50
```

### Passo 3: Obter Informações Taxonômicas

Consulte [NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy) para:
- `virus_tax_id`: ID taxonômico do vírus
- `virus_species_tax_id`: ID taxonômico da espécie

### Passo 4: Testar

```bash
vqc get-nextclade-datasets --cores 2
```

---

## Adicionando Dataset do GitHub

Para datasets personalizados hospedados no GitHub:

### Passo 1: Preparar Repositório

Seu repositório deve conter:

```
seu-repositorio/
└── caminho/do/dataset/
    ├── reference.fasta
    ├── genome_annotation.gff3
    ├── tree.json (opcional)
    ├── pathogen.json
    └── sequences.fasta (opcional)
```

### Passo 2: Adicionar ao datasets.yml

```yaml
github:
  meu-virus-customizado:
    repository: "seu-usuario/seu-repositorio"
    dataset: "caminho/dentro/do/repo"
    tag: "main"
    virus_name: "Meu Vírus Customizado"
    virus_tax_id: 123456
    virus_species: "Nome da Espécie"
    virus_species_tax_id: 789012
    segment: "Unsegmented"
    ncbi_id: "NC_XXXXXX.X"
    target_gene: "gene1"
    target_regions: ["gene1", "gene2"]
    private_mutation_total_threshold: 40
```

### Passo 3: Testar

```bash
vqc get-nextclade-datasets --cores 2
```

Verifique se o dataset foi baixado em `datasets/external_datasets/meu-virus-customizado/`.
