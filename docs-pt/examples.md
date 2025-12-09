# Exemplos Práticos

## Exemplo 1: Análise Básica de Dengue

```bash
# 1. Preparar ambiente (executar uma vez)
vqc get-nextclade-datasets
vqc get-blast-database

# 2. Executar análise
vqc run-from-fasta \
  --sequences-fasta amostras_dengue.fasta \
  --output-dir qc_dengue \
  --output-file resultados_dengue.tsv \
  --cores 4
```

## Exemplo 2: Influenza (Múltiplos Segmentos)

```bash
vqc run-from-fasta \
  --sequences-fasta amostras_influenza.fasta \
  --output-dir qc_influenza \
  --output-file influenza_qc.tsv \
  --cores 8
```

O ViralQC identifica automaticamente os diferentes segmentos (HA, NA, PB1, etc.).

## Exemplo 3: Análise Metagenômica

Para vírus desconhecidos, use parâmetros relaxados:

```bash
vqc run-from-fasta \
  --sequences-fasta metagenoma.fasta \
  --output-dir qc_metagenoma \
  --blast-pident 30 \
  --blast-qcov 30 \
  --blast-task dc-megablast \
  --ns-min-score 0.05
```

## Exemplo 4: Banco de Dados Reproduzível

Crie um banco BLAST com data específica de release:

```bash
# Criar banco com sequências até uma data específica
vqc get-blast-database --release-date 2023-06-15

# Executar análise com este banco
vqc run-from-fasta \
  --sequences-fasta amostras.fasta \
  --blast-database datasets/blast.fasta
```

## Exemplo 5: Saída em JSON

```bash
vqc run-from-fasta \
  --sequences-fasta amostras.fasta \
  --output-file resultados.json \
  --cores 4
```

## Exemplo 6: Busca BLAST Sensível

Para vírus divergentes:

```bash
vqc run-from-fasta \
  --sequences-fasta amostras.fasta \
  --blast-task blastn \
  --blast-pident 70 \
  --blast-evalue 1e-5 \
  --cores 4
```
