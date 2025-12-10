# Troubleshooting

## Problem 1: `vqc` command not found

**Symptom:**
```
bash: vqc: command not found
```

**Solution:**
```bash
# Activate environment
micromamba activate viralQC

# Verify installation
pip show viralQC

# Reinstall if needed
pip install --force-reinstall viralQC
```

## Problem 2: Error downloading datasets

**Symptom:**
```
ERROR: Failed to retrieve Nextclade public datasets
```

**Solutions:**
1. Check internet connection
2. Verify dataset path in `datasets.yml`
3. Update tag to most recent version

## Problem 3: BLAST database not found

**Symptom:**
```
ERROR: BLAST database not found at datasets/blast.fasta
```

**Solution:**
```bash
vqc get-blast-database
```

## Problem 4: No sequences mapped

**Symptom:** All sequences show as "Unclassified"

**Solutions:**
1. Check input sequence quality
2. Use relaxed BLAST parameters:
   ```bash
   --blast-pident 70 --blast-qcov 50
   ```
3. Add missing datasets to `datasets.yml`

## Problem 5: Memory errors

**Symptom:**
```
Killed/Core dump
```

**Solutions:**
1. Reduce cores: `--cores 2`
2. Split FASTA into smaller files
3. Use machine with more RAM

## Problem 6: Permission denied

**Solution:**
```bash
vqc run-from-fasta \
  --sequences-fasta seqs.fasta \
  --output-dir ~/my_results
```

## Problem 7: GitHub datasets not downloaded

**Solutions:**
1. Check repository format: `user/repository`
2. Verify tag/branch exists
3. Check GitHub connectivity: `ping github.com`

## Getting Help

- [GitHub Issues](https://github.com/InstitutoTodosPelaSaude/viralQC/issues)
- Check logs in `output/logs/` directory
