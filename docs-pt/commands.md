# Comandos e Uso

O ViralQC fornece três comandos principais através da CLI (`vqc`).

## get-nextclade-datasets

Baixa e configura os datasets do Nextclade localmente.

```{important}
Este comando deve ser executado **pelo menos uma vez** antes de usar o `run`.
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
| `--verbose` | Booleano | `False` | Exibir logs do snakemake |

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
Este comando deve ser executado **pelo menos uma vez** antes de usar o `run`.
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
| `--verbose` | Booleano | `False` | Exibir logs do snakemake |

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

## run

Comando principal de análise. Identifica vírus, realiza controle de qualidade e extrai regiões-alvo.

### Uso

```bash
vqc run --input minhas_sequencias.fasta
```

### Parâmetros Obrigatórios

| Parâmetro | Tipo | Descrição |
|-----------|------|-----------|
| `--input` | String | Caminho para o arquivo FASTA de entrada |

### Parâmetros de Saída

| Parâmetro | Tipo | Padrão | Descrição |
|-----------|------|--------|-----------|
| `--output-dir` | String | `outputs` | Diretório de trabalho. Os resultados serão armazenados em um subdiretório `outputs/` dentro desta pasta. |
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
vqc run --input seqs.fasta

# Busca mais sensível para vírus distantes
vqc run --input seqs.fasta --blast-task dc-megablast

# BLASTN tradicional para sequências divergentes
vqc run --input seqs.fasta --blast-task blastn
```

### Parâmetros de Sistema

| Parâmetro | Tipo | Padrão | Descrição |
|-----------|------|--------|-----------|
| `--cores` | Integer | `1` | Número de threads/cores |
| `--verbose` | Booleano | `False` | Exibir logs do snakemake |

### Exemplo Completo

```bash
vqc run \
  --input amostras.fasta \
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


## API

Como o propósito do viralQC é sua integração com banco de dados genômicos de vírus, é possível integrar o módulo de análise no código de outras aplicações.

Isso pode ser feito importando a classe `RunAnalysis` do módulo `viralqc.core.run_analysis`. Essa classe possui o método `run` que executa a análise de qualidade de um genoma viral, recebendo como parâmetro o caminho para o arquivo FASTA contendo as sequências a serem analisadas. Outros parametros podem ser informados de forma otimizada.

### Uso

```python
from viralqc.core.run_analysis import RunAnalysis

input_file = "seqs.fasta"
output_directory = "results"
output_file = "results.json"
run_analysis = RunAnalysis()

snakemake_response = run_analysis.run(
        sequences_fasta=input_file,
        output_dir=output_directory,
        output_file=output_file
)
```

Ou uma abordagem flexível

```python
from viralqc.core.run_analysis import RunAnalysis

input_file = "seqs.fasta"
output_directory = "results"
output_file = "results.json"
run_analysis = RunAnalysis()

snakemake_response = run_analysis.run(
        sequences_fasta=input_file,
        output_dir=output_directory,
        output_file=output_file,
        cores=2,
        datasets_local_path="datasets",
        nextclade_sort_min_score=0.1,
        nextclade_sort_min_hits=10,
        blast_database="datasets/blast.fasta",
        blast_database_metadata="datasets/blast.tsv",
        blast_identity_threshold=0,
        blast_evalue=0.01,
        blast_qcov=0,
        blast_task="blastn"
)
```

Para checar os resultados:

```python
if snakemake_response.status == 200:
    results_data = snakemake_response.get_results()
    for seq_result in results_data:
        virus = seq_result.get("virus")
        quality = seq_result.get("genomeQuality")
        coverage = seq_result.get("coverage")
        print(virus, quality, coverage)
else:
    raise Exception(snakemake_response.format_log())
```

### Atributos e Métodos

O método `run` retorna um objeto `SnakemakeResponse` que possui os atributos:

| Atributo | Tipo | Descrição |
|----------|------|-----------|
| run_id | str | ID da execução |
| status | RunStatus | Status da execução, que pode ser 200 (sucesso) ou 500 (falha) |
| log_path | str | Caminho para o arquivo de log |
| results_path | str | Caminho para o arquivo de resultados |

E os seguintes métodos:

| Método | Descrição |
|--------|-----------|
| `format_log()` | Retorna o conteúdo do arquivo de log formatado |
| `get_results()` | Retorna o conteúdo do arquivo de resultados em formato de dicionário |

