# Comandos e Uso

O ViralQC fornece três comandos principais através da CLI (`vqc`).

## get-nextclade-datasets

Baixa e configura os datasets do Nextclade localmente.

```{important}
Este comando deve ser executado **pelo menos uma vez** antes de usar o `run-from-fasta`.
```

### Uso

```bash
vqc get-nextclade-datasets --cores 2
```

### Parâmetros

| Parâmetro | Tipo | Padrão | Descrição |
|-----------|------|--------|-----------|
| `--datasets-dir` | String | `datasets` | Diretório para armazenar datasets |
| `--cores` | Integer | `1` | Número de threads/cores |

### Estrutura de Saída

```
datasets/
├── nextclade_data/
│   ├── denv1/
│   ├── denv2/
│   └── ...
├── external_datasets/
│   └── zikav/
└── external_datasets_minimizers.json
```

---

## get-blast-database

Cria um banco de dados BLAST local contendo genomas virais do NCBI RefSeq.

```{important}
Este comando deve ser executado **pelo menos uma vez** antes de usar o `run-from-fasta`.
```

### Uso

```bash
vqc get-blast-database
```

### Parâmetros

| Parâmetro | Tipo | Padrão | Descrição |
|-----------|------|--------|-----------|
| `--output-dir` | String | `datasets` | Diretório para o banco BLAST |
| `--release-date` | String | `None` | Filtrar sequências por data de release (YYYY-MM-DD). Apenas sequências liberadas até esta data serão incluídas |
| `--cores` | Integer | `1` | Número de threads/cores |

### Filtragem por Data de Release

O parâmetro `--release-date` permite criar um banco de dados BLAST reproduzível filtrando sequências pela data de release do NCBI:

```bash
# Criar banco com sequências liberadas até 15 de junho de 2023
vqc get-blast-database --release-date 2023-06-15
```

**Comportamento:**
- Quando `--release-date` é fornecido:
  - Apenas sequências com `release_date <= data_especificada` são incluídas
  - A data especificada é usada como identificador de versão do banco
- Quando não fornecido:
  - Todas as sequências RefSeq disponíveis são incluídas
  - A data atual é usada como identificador de versão

Isso é útil para:
- **Reprodutibilidade**: Recriar o mesmo banco em diferentes momentos
- **Auditoria**: Rastrear quais sequências estavam disponíveis em uma data específica
- **Estudos comparativos**: Analisar como os resultados mudam com atualizações do banco

### Versão do Banco

A versão do banco é registrada no arquivo de metadados `blast.tsv`:
- Formato: `ncbi-refseq-virus_YYYY-MM-DD`
- Usa o valor de `--release-date` se fornecido, caso contrário a data atual

### Estrutura de Saída

```
datasets/
├── blast.fasta          # Sequências de referência
├── blast.fasta.ndb      # Arquivos do banco BLAST
├── blast.fasta.nhr
├── blast.fasta.nin
├── blast.fasta.nsq
├── blast.tsv            # Metadados com versão
└── blast_gff/           # Arquivos GFF3
```

---

## run-from-fasta

Comando principal de análise. Identifica vírus, realiza controle de qualidade e extrai regiões-alvo.

### Uso

```bash
vqc run-from-fasta --sequences-fasta minhas_sequencias.fasta
```

### Parâmetros Obrigatórios

| Parâmetro | Tipo | Descrição |
|-----------|------|-----------|
| `--sequences-fasta` | String | Caminho para o arquivo FASTA de entrada |

### Parâmetros de Saída

| Parâmetro | Tipo | Padrão | Descrição |
|-----------|------|--------|-----------|
| `--output-dir` | String | `output` | Diretório para arquivos de saída |
| `--output-file` | String | `results.tsv` | Arquivo de resultados (`.tsv`, `.csv`, `.json`) |

### Parâmetros de Datasets

| Parâmetro | Tipo | Padrão | Descrição |
|-----------|------|--------|-----------|
| `--datasets-dir` | String | `datasets` | Caminho para diretório de datasets |

### Parâmetros do Nextclade Sort

| Parâmetro | Tipo | Padrão | Descrição |
|-----------|------|--------|-----------|
| `--ns-min-score` | Float | `0.1` | Score mínimo para match válido |
| `--ns-min-hits` | Integer | `10` | Hits mínimos para match válido |

### Parâmetros do BLAST

| Parâmetro | Tipo | Padrão | Descrição |
|-----------|------|--------|-----------|
| `--blast-database` | String | `datasets/blast.fasta` | Caminho para o banco BLAST |
| `--blast-database-metadata` | String | `datasets/blast.tsv` | Caminho para metadados |
| `--blast-pident` | Integer | `80` | Identidade mínima (0-100) |
| `--blast-evalue` | Float | `1e-10` | E-value máximo |
| `--blast-qcov` | Integer | `80` | Cobertura mínima da query (0-100) |
| `--blast-task` | String | `megablast` | Tipo de tarefa BLAST |

### Tipos de Tarefa BLAST

O parâmetro `--blast-task` controla a sensibilidade do algoritmo:

| Tarefa | Descrição | Caso de Uso |
|--------|-----------|-------------|
| `megablast` | Sequências muito similares (padrão) | Rápido, mesma espécie |
| `dc-megablast` | Discontiguous megablast | Entre espécies, mais sensível |
| `blastn` | BLASTN tradicional | Sequências distantes |
| `blastn-short` | Sequências curtas | Sequências < 50 bp |

**Exemplos:**

```bash
# Padrão (megablast) - rápido, para sequências similares
vqc run-from-fasta --sequences-fasta seqs.fasta

# Busca mais sensível para vírus distantes
vqc run-from-fasta --sequences-fasta seqs.fasta --blast-task dc-megablast

# BLASTN tradicional para sequências divergentes
vqc run-from-fasta --sequences-fasta seqs.fasta --blast-task blastn
```

### Parâmetros de Sistema

| Parâmetro | Tipo | Padrão | Descrição |
|-----------|------|--------|-----------|
| `--cores` | Integer | `1` | Número de threads/cores |

### Exemplo Completo

```bash
vqc run-from-fasta \
  --sequences-fasta amostras.fasta \
  --output-dir resultados \
  --output-file relatorio.tsv \
  --blast-pident 75 \
  --blast-task dc-megablast \
  --cores 8
```

### Fluxo de Análise

1. **Nextclade Sort**: Mapeia sequências para datasets locais
2. **Análise BLAST**: Identifica sequências não mapeadas
3. **Nextclade Run**: Análise de controle de qualidade
4. **Pós-processamento**: Combina e pontua resultados
5. **Extração de Regiões**: Extrai regiões-alvo baseado em qualidade
