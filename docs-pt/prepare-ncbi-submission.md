# Preparando Submissões ao NCBI

O `viralQC` inclui o comando `prepare-ncbi-submission` para organizar sequências e gerar CSVs de metadados formatados para submissão ao NCBI.

O comando contém dois sub-comandos:
1. `virus`: Agrupa sequências por vírus (ou tipos/segmentos).
2. `sample`: Agrupa sequências por vírus, mas apenas para IDs de amostras individuais.

Ambos os comandos requerem os arquivos de saída gerados pelo comando `viralQC run`:
* `--results`: O arquivo principal de resultados do ViralQC (ex: `results.tsv`).
* `--sequences-vqc`: O arquivo FASTA de saída com as regiões alvo (ex: `sequences_target_regions.fasta`) gerado pelo `viralQC run`. Este arquivo contém sequências alvo processadas pelo ViralQC e tem prioridade.
* `--sequences-input`: O arquivo FASTA de entrada original que você passou para o `viralQC run`.
  *(Nota: O comando prioriza sequências encontradas em `--sequences-vqc`. Se uma sequência foi filtrada/descartada pelo VQC (por falta de informação de qualidade), mas você ainda deseja submetê-la, ela será extraída do `--sequences-input`.)*
* `--output-prefix` (opcional, padrão `ncbi_submission`): Prefixo usado para os diretórios de saída gerados.

## Agrupamento por Vírus (`virus`)

Ao usar o sub-comando `virus`, as sequências são agrupadas em diretórios específicos para o vírus ou segmento.

```bash
vqc prepare-ncbi-submission virus [SUBCOMMAND] --results results.tsv --sequences sequences_target_regions.fasta
```

### Vírus Suportados

Para Dengue, Influenza, Norovírus e SARS-CoV-2, o NCBI possui requisitos específicos de submissão. Cada um exige colunas específicas no arquivo `metadata.csv`, e estes sub-comandos formatam o CSV de metadados adequadamente.

#### SARS-CoV-2
```bash
vqc prepare-ncbi-submission virus sars-cov-2 --results results.tsv --sequences-vqc sequences_target_regions.fasta --sequences-input original_sequences.fasta --metadata input_metadata.csv
```

Organiza todas as sequências de SARS-CoV-2 em `ncbi_submission_SARS-CoV-2/`.

#### Dengue
```bash
vqc prepare-ncbi-submission virus dengue --results results.tsv --sequences-vqc sequences_target_regions.fasta --sequences-input original_sequences.fasta --metadata input_metadata.csv
```

Organiza sequências por tipo, criando diretórios como `ncbi_submission_Dengue1/`, `ncbi_submission_Dengue2/`, etc., dependendo dos sorotipos identificados na análise do ViralQC.

#### Influenza
```bash
vqc prepare-ncbi-submission virus influenza --results results.tsv --sequences-vqc sequences_target_regions.fasta --sequences-input original_sequences.fasta --metadata input_metadata.csv
```

Organiza sequências por tipo com subdiretórios para cada segmento, criando diretórios como `ncbi_submission_InfluenzaA/HA/`, `ncbi_submission_InfluenzaA/NA/`, `ncbi_submission_InfluenzaB/HA/`, etc.

#### Norovírus
```bash
vqc prepare-ncbi-submission virus norovirus --results results.tsv --sequences-vqc sequences_target_regions.fasta --sequences-input original_sequences.fasta --metadata input_metadata.csv
```

Organiza sequências por genogrupo, criando subdiretórios como `ncbi_submission_Norovirus/GI/`, `ncbi_submission_Norovirus/GII/`, etc., dependendo dos genogrupos identificados na análise do ViralQC. Os genogrupos suportados vão de GI a GVI.

#### Vírus Personalizados (Custom Viruses)
```bash
vqc prepare-ncbi-submission virus custom \
    --virus-name "Respiratory syncytial virus A" \
    --results results.tsv \
    --sequences-vqc sequences_target_regions.fasta \
    --sequences-input original_sequences.fasta \
    --metadata input_metadata.csv
```

Organiza todas as sequências que correspondem ao `--virus-name` fornecido em um único diretório. Vírus não padronizados recebem automaticamente o qualificador `[Organism=...]` em seus cabeçalhos FASTA com base na `virus_species` identificada. O nome fornecido deve ser o mesmo presente no campo `virus` do arquivo de resultados.

Se o vírus possui segmentos anotados (ex: S, M, L), você pode passar a flag `--split-by-segments` para organizar as sequências em subdiretórios por segmento (ex: `ncbi_submission_Oropouche_virus/S/`, `ncbi_submission_Oropouche_virus/M/`, `ncbi_submission_Oropouche_virus/L/`). Cada subdiretório conterá seus próprios arquivos `sequences.fasta`, `metadata.tsv`, `annotation.tbl` e arquivos de log. Sem esta flag, todas as sequências são colocadas em um único diretório.

Além disso, para vírus personalizados, você pode passar `--tbl-dir` apontando para uma pasta com arquivos de anotação `.tbl` por amostra. Estes serão concatenados em um único arquivo `annotation.tbl` junto com as sequências FASTA. Se você não fornecer esta opção, o comando tentará encontrar o `.tbl` com base no campo `tbl_path` no arquivo de resultados.

## Agrupamento por Amostra (`sample`)

Se você preferir organizar as submissões apenas para amostras específicas, ou para todas as amostras independentemente do vírus, use o sub-comando `sample`.

```bash
vqc prepare-ncbi-submission sample \
    --sample <sample_id_1> \
    --sample <sample_id_2> \
    --results results.tsv \
    --sequences-vqc sequences_target_regions.fasta \
    --sequences-input original_sequences.fasta \
    --metadata input_metadata.csv
```
Você pode passar múltiplas opções `--sample`, ou fornecer um arquivo de texto com um ID por linha via `--sample-ids samples.txt`.

Para processar **todas** as amostras presentes no arquivo de resultados, basta usar:
```bash
vqc prepare-ncbi-submission sample --sample all ...
```

Para processar amostras baseadas em uma lista de IDs de amostras, use a opção `--sample-ids`:
```bash
vqc prepare-ncbi-submission sample --sample-ids samples.txt ...
```

Isso cria um diretório por vírus (ex: `ncbi_submission_Dengue1/`). Se uma sequência foi filtrada ou faltaram dados (comumente arquivos tbl), isso é informado em um arquivo `[prefix]_skipped.tsv`.

## CSV de Metadados

Para todos os comandos de submissão, você pode (ou deve, para vírus predefinidos) fornecer um arquivo CSV de metadados de entrada via `--metadata`.

### Colunas de Entrada

O CSV de entrada pode conter as seguintes colunas. Sua necessidade depende do tipo de vírus sendo processado:

| Coluna | Descrição | Dengue, Influenza e Norovírus | SARS-CoV-2 | Vírus Personalizados |
|--------|-----------|-------------------------------|------------|----------------------|
| `Sequence_ID` | Deve corresponder exatamente ao `seqName` como aparece no arquivo de resultados e nos cabeçalhos FASTA. **Deve ter menos de 25 caracteres.** | **Obrigatório** | **Obrigatório** | **Obrigatório** |
| `geo_loc_name` | A localização geográfica da amostra (ex: País). | **Obrigatório** | **Obrigatório** | Opcional |
| `host` | O hospedeiro natural do vírus (ex: `Homo sapiens`). Não use caracteres especiais. | **Obrigatório** | **Obrigatório** | Opcional |
| `isolate` | O nome do isolado ou string identificadora. | **Obrigatório** | **Obrigatório** | Opcional |
| `collection-date` | A data em que a amostra foi coletada, normalmente usando o formato `AAAA-MM-DD`. | **Obrigatório** | **Obrigatório** | Opcional |
| `isolation-source` | O material de origem da amostra (ex: `Serum`, `Swab`). | **Obrigatório** | *Ignorado* | Opcional |

*Nota: Para vírus personalizados, a criação do arquivo `--metadata` em si é completamente opcional. Se fornecido, apenas `Sequence_ID` deve estar presente, e as outras colunas de dados serão incluídas se você as adicionar.*

### Formato de Saída

O comando `prepare-ncbi-submission` pegará seu CSV de entrada e gerará um arquivo de metadados final dentro de cada diretório de submissão. Para vírus predefinidos (SARS-CoV-2, Dengue, Influenza, Norovírus), um `metadata.csv` delimitado por vírgula é gerado. Para vírus personalizados, um `metadata.tsv` delimitado por tabulação é gerado. A ferramenta enriquece automaticamente este arquivo com dados taxonômicos e de tipagem derivados dos resultados do `viralQC`:

* **Influenza**: Adiciona `serotype` (ex: `H1N1` ou `H3N2`) extraído diretamente da classificação.
* **Dengue**: Adiciona `genotype` (ex: `1`, `2`, `3` ou `4`) e `serotype` (atribuição detalhada de clado).
* **Norovírus**: Adiciona `genotype` (ex: `GII`, `GII.17`) extraído da classificação viral. Nenhuma coluna `serotype` é incluída.
* **SARS-CoV-2**: Remove `isolation-source` e `serotype`, pois normalmente não são incluídos em submissões de SC2 ao NCBI.
* **Vírus Personalizados**: Renomeia as colunas para nomes compatíveis com o NCBI: `geo_loc_name` → `Country (geo_loc_name)`, `host` → `Host`, `isolate` → `Isolate`, `collection-date` → `Collection_date`, `isolation-source` → `Isolation_source`. Nenhuma coluna `Organism` é adicionada (a informação do organismo é incluída nos cabeçalhos FASTA como `[Organism=...]`).

## Cabeçalhos FASTA e Anotações

Os cabeçalhos FASTA são gerenciados cuidadosamente durante a organização:
* Sequências que não atingem os limites de qualidade (ex: comprimento < 150nt, ou conteúdo de N ≥ 50%) são omitidas do FASTA e registradas em `dropped_sequences.tsv`.
* Caracteres não seguros em nomes de sequências (não-ASCII ou pipes) são sanitizados para sublinhados para compatibilidade com o NCBI, com as traduções registradas em `renamed_headers.tsv`.
* Espaços e colchetes são preservados corretamente, permitindo que qualificadores de recursos padrão do NCBI, como `[Organism=...]`, funcionem como pretendido para vírus não padronizados.

## Divisão em Lotes (Batch Splitting)

O NCBI limita submissões a 3.000 sequências por arquivo. Quando um grupo de vírus excede esse limite, os arquivos `sequences.fasta` e de metadados são automaticamente divididos em lotes numerados:

* `sequences.1.fasta`, `metadata.1.csv` (ou `metadata.1.tsv` para vírus personalizados) — primeiras 2.999 sequências
* `sequences.2.fasta`, `metadata.2.csv` (ou `metadata.2.tsv`) — próximas 2.999 sequências
* …e assim por diante.

Se um grupo tiver 2.999 ou menos sequências, os arquivos são escritos normalmente sem nenhum sufixo.
