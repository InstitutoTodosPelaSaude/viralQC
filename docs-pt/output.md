# Estrutura de Saída

Quando você executa `run-from-fasta`, o ViralQC cria a seguinte estrutura:

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

## Arquivo Principal de Resultados

O `results.tsv` contém:

### Colunas de Identificação

| Coluna | Descrição |
|--------|-----------|
| `seqName` | Nome da sequência do FASTA |
| `virus` | Nome do vírus identificado |
| `virus_tax_id` | ID taxonômico do NCBI |
| `segment` | Segmento genômico |
| `clade` | Clado filogenético |
| `datasetVersion` | Versão do dataset usado |

### Colunas de Qualidade (ViralQC)

| Coluna | Descrição |
|--------|-----------|
| `genomeQuality` | Qualidade geral (A, B, C, D) |
| `genomeQualityScore` | Score numérico (0-24) |
| `missingDataQuality` | Qualidade de dados ausentes |
| `privateMutationsQuality` | Qualidade de mutações privadas |
| `cdsCoverageQuality` | Qualidade de cobertura de CDS |
| `targetRegionsQuality` | Qualidade das regiões-alvo |

### Colunas de Mutações

| Coluna | Descrição |
|--------|-----------|
| `totalSubstitutions` | Total de substituições nucleotídicas |
| `totalDeletions` | Total de deleções |
| `totalInsertions` | Total de inserções |
| `coverage` | Cobertura do genoma (0.0-1.0) |

## Arquivos de Regiões-Alvo

### sequences_target_regions.bed

```
seq1    94      2419    C,prM,E
seq2    0       10735   genome
```

### sequences_target_regions.fasta

Sequências extraídas das regiões que atendem aos critérios de qualidade.
