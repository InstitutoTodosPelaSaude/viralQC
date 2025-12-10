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

### Arquivo Principal: results.tsv (ou .csv, .json)

Este é o arquivo que contém os resultados consolidados de todas as análises:

##### 1. Identificação da Sequência

| Coluna | Tipo | Descrição |
|--------|------|-----------|
| `seqName` | String | Nome da sequência no arquivo FASTA de entrada |
| `virus` | String | Nome do vírus identificado |
| `virus_tax_id` | Integer | ID taxonômico do vírus no NCBI Taxonomy |
| `virus_species` | String | Nome da espécie viral |
| `virus_species_tax_id` | Integer | ID taxonômico da espécie |
| `segment` | String | Segmento genômico (ex: "HA", "NA", "Unsegmented") |
| `ncbi_id` | String | Accession do genoma de referência no NCBI |
| `dataset` | String | Identificador do dataset usado |
| `datasetVersion` | String | Versão/tag do dataset |
| `clade` | String | Clado filogenético (quando disponível) |

##### 2. Métricas de Qualidade (Sistema A-D do ViralQC)

| Coluna | Tipo | Descrição |
|--------|------|-----------|
| `genomeQuality` | String | Qualidade geral do genoma (A, B, C ou D) |
| `genomeQualityScore` | Float | Score numérico normalizado (0-24) |
| `missingDataQuality` | String | Qualidade baseada em dados ausentes (A, B, C ou D) |
| `privateMutationsQuality` | String | Qualidade baseada em mutações privadas (A, B, C ou D) |
| `mixedSitesQuality` | String | Qualidade baseada em sítios mistos (A, B, C ou D) |
| `snpClustersQuality` | String | Qualidade baseada em clusters de SNPs (A, B, C ou D) |
| `frameShiftsQuality` | String | Qualidade baseada em frameshifts (A, B, C ou D) |
| `stopCodonsQuality` | String | Qualidade baseada em códons de parada (A, B, C ou D) |
| `cdsCoverageQuality` | String | Qualidade de cobertura por CDS (ex: "E: A, prM: B") |
| `targetRegionsQuality` | String | Qualidade das regiões-alvo (A, B, C, D ou vazio) |
| `targetGeneQuality` | String | Qualidade do gene-alvo (A, B, C, D ou vazio) |

##### 3. Métricas de Qualidade (Nextclade)

| Coluna | Tipo | Descrição |
|--------|------|-----------|
| `qc.overallScore` | Float | Score geral de qualidade do Nextclade |
| `qc.overallStatus` | String | Status de qualidade do Nextclade (good, mediocre, bad) |
| `qc.privateMutations.total` | Integer | Total de mutações privadas (Nextclade) |
| `qc.privateMutations.score` | Float | Score de mutações privadas (Nextclade) |
| `qc.privateMutations.status` | String | Status de mutações privadas (Nextclade) |
| `qc.missingData.score` | Float | Score de dados ausentes (Nextclade) |
| `qc.missingData.status` | String | Status de dados ausentes (Nextclade) |
| `qc.mixedSites.totalMixedSites` | Integer | Total de sítios mistos (Nextclade) |
| `qc.mixedSites.score` | Float | Score de sítios mistos (Nextclade) |
| `qc.mixedSites.status` | String | Status de sítios mistos (Nextclade) |
| `qc.snpClusters.totalSNPs` | Integer | Total de SNPs em clusters (Nextclade) |
| `qc.snpClusters.score` | Float | Score de clusters de SNPs (Nextclade) |
| `qc.snpClusters.status` | String | Status de clusters de SNPs (Nextclade) |
| `qc.frameShifts.totalFrameShifts` | Integer | Total de frameshifts (Nextclade) |
| `qc.frameShifts.score` | Float | Score de frameshifts (Nextclade) |
| `qc.frameShifts.status` | String | Status de frameshifts (Nextclade) |
| `qc.stopCodons.totalStopCodons` | Integer | Total de códons de parada (Nextclade) |
| `qc.stopCodons.score` | Float | Score de códons de parada (Nextclade) |
| `qc.stopCodons.status` | String | Status de códons de parada (Nextclade) |

##### 4. Cobertura e Regiões

| Coluna | Tipo | Descrição |
|--------|------|-----------|
| `coverage` | Float | Cobertura do genoma (0.0 a 1.0) |
| `cdsCoverage` | String | Cobertura de cada CDS (formato: "gene1: 0.98, gene2: 1.0") |
| `targetRegionsCoverage` | String | Cobertura das regiões-alvo definidas em `target_regions` |
| `targetGeneCoverage` | String | Cobertura do gene-alvo definido em `target_gene` |
| `targetRegions` | String | Lista de regiões-alvo (separadas por \|) |
| `targetGene` | String | Nome do gene-alvo principal |

##### 5. Mutações Nucleotídicas

| Coluna | Tipo | Descrição |
|--------|------|-----------|
| `totalSubstitutions` | Integer | Total de substituições nucleotídicas |
| `totalDeletions` | Integer | Total de deleções nucleotídicas |
| `totalInsertions` | Integer | Total de inserções nucleotídicas |
| `totalFrameShifts` | Integer | Total de frameshift mutations |
| `totalMissing` | Integer | Total de nucleotídeos faltantes (N's ou gaps) |
| `totalNonACGTNs` | Integer | Total de caracteres não-ACGTN |
| `substitutions` | String | Lista de substituições (formato: gene:pos:ref>alt) |
| `deletions` | String | Lista de deleções |
| `insertions` | String | Lista de inserções |
| `frameShifts` | String | Lista de frameshifts |
| `alignmentScore` | Float | Score do alinhamento |

##### 6. Mutações de Aminoácidos

| Coluna | Tipo | Descrição |
|--------|------|-----------|
| `totalAminoacidSubstitutions` | Integer | Total de substituições de aminoácidos |
| `totalAminoacidDeletions` | Integer | Total de deleções de aminoácidos |
| `totalAminoacidInsertions` | Integer | Total de inserções de aminoácidos |
| `totalUnknownAa` | Integer | Total de aminoácidos desconhecidos |
| `aaSubstitutions` | String | Lista de substituições de aminoácidos |
| `aaDeletions` | String | Lista de deleções de aminoácidos |
| `aaInsertions` | String | Lista de inserções de aminoácidos |

##### 7. Mutações Privadas (Detalhadas)

| Coluna | Tipo | Descrição |
|--------|------|-----------|
| `privateNucMutations.totalPrivateSubstitutions` | Integer | Total de substituições privadas |
| `privateNucMutations.totalLabeledSubstitutions` | Integer | Total de mutações privadas conhecidas/catalogadas |
| `privateNucMutations.totalUnlabeledSubstitutions` | Integer | Total de mutações privadas não catalogadas |
| `privateNucMutations.totalReversionSubstitutions` | Integer | Total de reversões (mutações que revertem à referência ancestral) |

**Nota sobre formatos de saída:**
- **TSV/CSV**: Todas as colunas são strings ou valores numéricos
- **JSON**: Colunas como `cdsCoverage`, `cdsCoverageQuality` e `targetRegionsCoverage` são formatadas como arrays de objetos para facilitar parsing programático


## Arquivos de Regiões-Alvo

### sequences_target_regions.bed

```
seq1    94      2419    C,prM,E
seq2    0       10735   genome
```

### sequences_target_regions.fasta

Sequências extraídas das regiões que atendem aos critérios de qualidade.
