# Solução de Problemas

## Problema 1: Comando `vqc` não encontrado

**Sintoma:**
```
bash: vqc: command not found
```

**Solução:**
```bash
# Ativar ambiente
micromamba activate viralQC

# Verificar instalação
pip show viralQC

# Reinstalar se necessário
pip install --force-reinstall viralQC
```

## Problema 2: Erro ao baixar datasets

**Sintoma:**
```
ERROR: Failed to retrieve Nextclade public datasets
```

**Soluções:**
1. Verificar conexão com internet
2. Verificar caminho do dataset em `datasets.yml`
3. Atualizar tag para versão mais recente

## Problema 3: Banco BLAST não encontrado

**Sintoma:**
```
ERROR: BLAST database not found at datasets/blast.fasta
```

**Solução:**
```bash
vqc get-blast-database
```

## Problema 4: Nenhuma sequência mapeada

**Sintoma:** Todas as sequências aparecem como "Unclassified"

**Soluções:**
1. Verificar qualidade das sequências de entrada
2. Usar parâmetros BLAST relaxados:
   ```bash
   --blast-pident 70 --blast-qcov 50
   ```
3. Adicionar datasets faltantes ao `datasets.yml`

## Problema 5: Erros de memória

**Sintoma:**
```
Killed/Core dump
```

**Soluções:**
1. Reduzir cores: `--cores 2`
2. Dividir FASTA em arquivos menores
3. Usar máquina com mais RAM

## Problema 6: Permissão negada

**Solução:**
```bash
vqc run-from-fasta \
  --sequences-fasta seqs.fasta \
  --output-dir ~/meus_resultados
```

## Problema 7: Datasets do GitHub não baixados

**Soluções:**
1. Verificar formato do repositório: `usuario/repositorio`
2. Verificar se tag/branch existe
3. Verificar conectividade: `ping github.com`

## Problema 8: BLAST rejeita uma sequência como dado implausível

**Sintoma:** a execução para na regra `blast` e `output/logs/blast.log` contém:

```
BLAST query error: CFastaReader: Near line 2, there's a line that doesn't look like plausible data, but it's not marked as defline or comment.
```

**Causa:** o FASTA de entrada contém um registro composto inteiramente por
caracteres de gap de alinhamento (`-` ou `.`), tipicamente um consenso vazio
exportado como linha de alinhamento. O `blastn` recusa esse registro e, como a
regra roda sob `set -euo pipefail`, a falha abortava toda a execução.

**Solução:** atualize para a **v1.2.0 ou posterior**. Registros sem nenhuma base
utilizável passam a ser excluídos do Nextclade e do BLAST, listados em
`excluded_sequences.tsv` e reportados em `results.tsv` com
`inputSequenceStatus` igual a `no usable sequence`. Não é necessário alterar o
arquivo de entrada.

Em versões anteriores, remova os gaps e descarte os registros vazios antes de
executar (`-l 0` coloca cada sequência em uma única linha, o que o filtro exige):

```bash
seqtk seq -A -l 0 entrada.fasta | awk '/^>/{h=$0; next} {gsub(/[-.]/,"",$0); if (length($0)>0) print h"\n"$0}' > limpo.fasta
```

## Obtendo Ajuda

- [GitHub Issues](https://github.com/InstitutoTodosPelaSaude/viralQC/issues)
- Verifique logs no diretório `output/logs/`
