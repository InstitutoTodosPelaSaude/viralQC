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

#### Colunas Obrigatórias

| Coluna | Descrição | Dengue, Influenza e Norovírus | SARS-CoV-2 | Vírus Personalizados |
|--------|-----------|-------------------------------|------------|----------------------|
| `Sequence_ID` | Deve corresponder exatamente ao `seqName` como aparece no arquivo de resultados e nos cabeçalhos FASTA. **Deve ter menos de 25 caracteres.** | **Obrigatório** | **Obrigatório** | **Obrigatório** |
| `geo_loc_name` | A localização geográfica da amostra (ex: País). | **Obrigatório** | **Obrigatório** | Opcional |
| `host` | O hospedeiro natural do vírus (ex: `Homo sapiens`). Não use caracteres especiais. | **Obrigatório** | **Obrigatório** | Opcional |
| `isolate` | O nome do isolado ou string identificadora. | **Obrigatório** | **Obrigatório** | Opcional |
| `collection-date` | A data em que a amostra foi coletada, normalmente usando o formato `AAAA-MM-DD`. | **Obrigatório** | **Obrigatório** | Opcional |
| `isolation-source` | O material de origem da amostra (ex: `Serum`, `Swab`). | **Obrigatório** | *Ignorado* | Opcional |

*Nota: Para vírus personalizados, a criação do arquivo `--metadata` em si é completamente opcional. Se fornecido, apenas `Sequence_ID` deve estar presente, e as outras colunas de dados serão incluídas se você as adicionar.*

#### Colunas Opcionais — Vírus Padrão (Dengue, Influenza, Norovírus, SARS-CoV-2)

Qualquer um dos seguintes modificadores INSDC pode ser adicionado ao CSV de metadados. Se presente, será automaticamente incluído no `metadata.csv` de saída:

| Coluna (cabeçalho CSV) | Descrição |
|------------------------|----------|
| `altitude` | Altitude em metros acima ou abaixo do nível do mar onde a amostra foi coletada. |
| `collected_by` | Nome da pessoa que coletou a amostra. |
| `culture_collection` | Código da instituição e ID da cultura (formato: `inst:col:id`). |
| `haplotype` | Halotipo do organismo. |
| `lab_host` | Hospedeiro laboratorial usado para propagar o organismo. |
| `lat_lon` | Latitude e longitude em graus decimais (ex: `15.77 S 47.93 W`). |
| `note` | Qualquer informação adicional em texto livre sobre a sequência. |
| `segment` | Nome do segmento viral ou de fago sequenciado. |
| `sex` | Sexo do organismo do qual a sequência foi obtida. |
| `specimen_voucher` | Identificador institucional para o espécime fonte. |
| `strain` | Cepa do organismo. |
| `tissue_type` | Tipo de tecido do qual a sequência foi obtida. |

#### Colunas Opcionais — Vírus Personalizados

Para vírus personalizados, os metadados são um arquivo delimitado por tabulação e os nomes das colunas seguem a convenção Title_Case do BankIt. As seguintes colunas podem ser adicionadas ao CSV de entrada usando os nomes em minúsculas padrão — elas serão automaticamente renomeadas na saída:

| Coluna de entrada | Coluna de saída | Descrição |
|---|---|---|
| `altitude` | `Altitude` | Altitude em metros. |
| `bio_material` | `Bio_material` | Identificador de material biológico. |
| `breed` | `Breed` | Raça nomeada (geralmente para mamíferos domesticados). |
| `cell_line` | `Cell_line` | Linhagem celular da qual a sequência foi obtida. |
| `cell_type` | `Cell_type` | Tipo de célula. |
| `clone` | `Clone` | Clone do qual a sequência foi obtida. |
| `collected_by` | `Collected_by` | Pessoa que coletou a amostra. |
| `culture_collection` | `Culture_collection` | Identificador de coleção de cultura. |
| `dev_stage` | `Dev_stage` | Estágio de desenvolvimento do organismo. |
| `ecotype` | `Ecotype` | Ecotipo nomeado. |
| `fwd_primer_name` | `Fwd_primer_name` | Nome do primer PCR direto. |
| `fwd_primer_seq` | `Fwd_primer_seq` | Sequência do primer PCR direto. |
| `genotype` | `Genotype` | Genótipo do organismo. |
| `haplogroup` | `Haplogroup` | Haplogrupo do organismo. |
| `haplotype` | `Haplotype` | Halotipo do organismo. |
| `lab_host` | `Lab_host` | Hospedeiro laboratorial usado para propagar o organismo. |
| `lat_lon` | `Lat_Lon` | Latitude e longitude em graus decimais. |
| `note` | `Note` | Informação adicional em texto livre. |
| `rev_primer_name` | `Rev_primer_name` | Nome do primer PCR reverso. |
| `rev_primer_seq` | `Rev_primer_seq` | Sequência do primer PCR reverso. |
| `segment` | `Segment` | Segmento viral ou de fago sequenciado. |
| `serotype` | `Serotype` | Variedade sorológica. |
| `serovar` | `Serovar` | Variedade sorológica (procarioto). |
| `sex` | `Sex` | Sexo do organismo. |
| `specimen_voucher` | `Specimen_voucher` | Identificador do espécime. |
| `strain` | `Strain` | Cepa do organismo. |
| `sub_species` | `Sub_species` | Subespécie. |
| `tissue_lib` | `Tissue_lib` | Biblioteca de tecido. |
| `tissue_type` | `Tissue_type` | Tipo de tecido. |
| `variety` | `Variety` | Variedade do organismo. |

> **Nota:** Colunas não listadas acima são ignoradas silenciosamente e não aparecerão no arquivo de metadados de saída.

### Formato de Saída

O comando `prepare-ncbi-submission` pegará seu CSV de entrada e gerará um arquivo de metadados final dentro de cada diretório de submissão. Para vírus predefinidos (SARS-CoV-2, Dengue, Influenza, Norovírus), um `metadata.csv` delimitado por vírgula é gerado. Para vírus personalizados, um `metadata.tsv` delimitado por tabulação é gerado. A ferramenta enriquece automaticamente este arquivo com dados taxonômicos e de tipagem derivados dos resultados do `viralQC`:

* **Influenza**: Adiciona `serotype` (ex: `H1N1` ou `H3N2`) extraído diretamente da classificação.
* **Dengue**: Adiciona `genotype` (ex: `1`, `2`, `3` ou `4`) e `serotype` (atribuição detalhada de clado).
* **Norovírus**: Adiciona `genotype` (ex: `GII`, `GII.17`) extraído da classificação viral. Nenhuma coluna `serotype` é incluída.
* **SARS-CoV-2**: Remove `isolation-source` e `serotype`, pois normalmente não são incluídos em submissões de SC2 ao NCBI.
* **Vírus Personalizados**: Renomeia as colunas para nomes compatíveis com o NCBI: `geo_loc_name` → `Country (geo_loc_name)`, `host` → `Host`, `isolate` → `Isolate`, `collection-date` → `Collection_date`, `isolation-source` → `Isolation_source`. Nenhuma coluna `Organism` é adicionada (a informação do organismo é incluída nos cabeçalhos FASTA como `[Organism=...]`).

## Cabeçalhos FASTA e Anotações

Os cabeçalhos FASTA são gerenciados cuidadosamente durante a organização:
* Sequências que não atingem os limites de qualidade (ex: comprimento < 150nt, ou conteúdo de N ≥ 50%) são excluídas do FASTA. O motivo de cada exclusão é registrado no arquivo `summary.txt` dentro do diretório de saída correspondente.
* **IDs de sequências duplicados de forma case-insensitive são detectados e removidos automaticamente.** O NCBI trata IDs de sequências como case-insensitive — por exemplo, `SEQ001` e `seq001` são considerados idênticos pelo NCBI e causariam um erro no upload. Quando um conflito é detectado, a **segunda** ocorrência é descartada e o motivo é registrado no `summary.txt`:
  ```
  dropped: 1 sequence(s)
    - seq001: Case-insensitive duplicate of 'SEQ001' (NCBI treats IDs as case-insensitive)
  ```
* Caracteres não seguros em nomes de sequências (não-ASCII ou pipes) são sanitizados para sublinhados para compatibilidade com o NCBI, com as traduções registradas em `renamed_headers.tsv`.
* Espaços e colchetes são preservados corretamente, permitindo que qualificadores de recursos padrão do NCBI, como `[Organism=...]`, funcionem como pretendido para vírus não padronizados.

## Divisão em Lotes (Batch Splitting)

O NCBI limita submissões a 3.000 sequências por arquivo. Quando um grupo de vírus excede esse limite, os arquivos `sequences.fasta` e de metadados são automaticamente divididos em lotes numerados:

* `sequences.1.fasta`, `metadata.1.csv` (ou `metadata.1.tsv` para vírus personalizados) — primeiras 2.999 sequências
* `sequences.2.fasta`, `metadata.2.csv` (ou `metadata.2.tsv`) — próximas 2.999 sequências
* …e assim por diante.

Se um grupo tiver 2.999 ou menos sequências, os arquivos são escritos normalmente sem nenhum sufixo.

## API Python

A lógica de preparo também está disponível como uma classe Python, `PrepareSubmission`, para uso em scripts ou pipelines de terceiros. Em vez de ler metadados de um arquivo CSV, a classe aceita uma lista de dicionários Python.

### Instalação

A classe está disponível após instalar o pacote `viralqc`. Nenhuma dependência adicional é necessária.

```python
from viralqc.core import PrepareSubmission
```

### Construtor

```python
PrepareSubmission(
    viralqc_results,      # Path – arquivo de resultados do ViralQC (.tsv, .csv ou .json)
    viralqc_target_seq,   # Path – sequences_target_regions.fasta gerado pelo vqc run
    viralqc_input_seq,    # Path – FASTA de entrada original passado ao vqc run
    samples_metadata,     # list[dict] – metadados das amostras (veja abaixo)
    output_prefix="ncbi_submission",  # str  – prefixo para os diretórios de saída
    split_by_segments=False,          # bool – dividir vírus personalizados por segmento
    tbl_dir=None,                     # Path|None – pasta com arquivos .tbl por amostra
)
```

Cada dicionário em `samples_metadata` utiliza as seguintes chaves:

| Chave | Descrição | Vírus padrão | Vírus personalizados |
|-------|-----------|--------------|----------------------|
| `sample_id` | Deve corresponder exatamente ao `seqName`. **Máximo 24 caracteres.** | **Obrigatório** | **Obrigatório** |
| `country` | Localização geográfica (mapeado para `geo_loc_name`). | **Obrigatório** | Opcional |
| `host` | Hospedeiro natural (ex: `Homo sapiens`). | **Obrigatório** | Opcional |
| `isolate` | Nome ou identificador do isolado. | **Obrigatório** | Opcional |
| `collection-date` | Data de coleta (`AAAA-MM-DD`). | **Obrigatório** | Opcional |
| `isolation-source` | Material de origem (ex: `Soro`). | **Obrigatório** | Opcional |

### Métodos

#### `run_virus(virus="all", virus_name=None)`

Prepara pacotes de submissão agrupados por tipo de vírus. Equivalente a `vqc prepare-ncbi-submission virus <subcomando>`.

- `virus`: `"all"` (padrão), `"sars-cov-2"`, `"dengue"`, `"influenza"`, `"norovirus"` ou `"custom"`.
- `virus_name`: obrigatório quando `virus="custom"`.

#### `run_sample(samples=["all"])`

Prepara pacotes para amostras específicas ou todas as amostras. Equivalente a `vqc prepare-ncbi-submission sample`.

- `samples`: lista de IDs de amostras, ou `["all"]` para processar todas.

### Valor de Retorno

Ambos os métodos retornam uma lista de dicionários, uma entrada por diretório de saída gerado:

```python
[
    {
        "SARS-CoV-2": {
            "sequences":  [Path("ncbi_submission_SARS-CoV-2/sequences.fasta")],
            "metadata":   [Path("ncbi_submission_SARS-CoV-2/metadata.csv")],
            "log":         Path("ncbi_submission_SARS-CoV-2/summary.txt"),
        }
    },
    {
        "Oropouche virus": {
            "sequences":  [Path("ncbi_submission_Oropouche_virus/sequences.fasta")],
            "metadata":   [Path("ncbi_submission_Oropouche_virus/metadata.tsv")],
            "log":         Path("ncbi_submission_Oropouche_virus/summary.txt"),
            "annotation": [Path("ncbi_submission_Oropouche_virus/annotation.tbl")],
        }
    },
]
```

A chave `"annotation"` está presente apenas para vírus personalizados que possuem arquivos TBL.

Para vírus organizados em subdiretórios (segmentos do Influenza, genogrupos do Norovírus, vírus personalizados com `split_by_segments=True`), cada subdiretório gera sua própria entrada. O rótulo usa o formato `"Tipo/Subgrupo"`, ex: `"InfluenzaA/HA"` ou `"Norovirus/GII"`.

### Exemplos

#### Processar todos os vírus encontrados no arquivo de resultados

```python
from pathlib import Path
from viralqc.core import PrepareSubmission

ps = PrepareSubmission(
    viralqc_results=Path("results.tsv"),
    viralqc_target_seq=Path("sequences_target_regions.fasta"),
    viralqc_input_seq=Path("sequences.fasta"),
    samples_metadata=[
        {
            "sample_id": "S001",
            "country": "Brasil",
            "host": "Homo sapiens",
            "isolate": "isolate/S001/2024",
            "collection-date": "2024-01-01",
            "isolation-source": "Soro",
        },
    ],
)

resultados = ps.run_virus()  # processa todos os grupos de vírus
for entrada in resultados:
    for rotulo, arquivos in entrada.items():
        print(f"{rotulo}:")
        for seq in arquivos["sequences"]:
            print(f"  sequências → {seq}")
        for meta in arquivos["metadata"]:
            print(f"  metadados  → {meta}")
        if "annotation" in arquivos:
            for ann in arquivos["annotation"]:
                print(f"  anotação   → {ann}")
```

#### Processar apenas amostras específicas

```python
resultados = ps.run_sample(samples=["S001"])
```

#### Processar vírus personalizado com divisão por segmentos

```python
ps = PrepareSubmission(
    viralqc_results=Path("results.tsv"),
    viralqc_target_seq=Path("sequences_target_regions.fasta"),
    viralqc_input_seq=Path("sequences.fasta"),
    samples_metadata=[...],
    split_by_segments=True,
)
resultados = ps.run_virus(virus="custom", virus_name="Oropouche virus")
# Gera ncbi_submission_Oropouche_virus/S/, /M/, /L/
```
