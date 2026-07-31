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

## Problem 8: BLAST rejects a sequence as implausible data

**Symptom:** the run stops at the `blast` rule and `output/logs/blast.log`
contains:

```
BLAST query error: CFastaReader: Near line 2, there's a line that doesn't look like plausible data, but it's not marked as defline or comment.
```

**Cause:** the input FASTA contains a record made up entirely of alignment gap
characters (`-` or `.`), typically an empty consensus exported as an alignment
row. `blastn` refuses such a record, and because the rule runs under
`set -euo pipefail` the failure aborted the whole run.

**Solution:** upgrade to **v1.2.0 or later**. Records with no usable base are
now excluded from Nextclade and BLAST, listed in `excluded_sequences.tsv`, and
reported in `results.tsv` with `inputSequenceStatus` set to
`no usable sequence`. No change to the input file is required.

On earlier versions, strip the gaps and drop the empty records before running
(`-l 0` puts each sequence on a single line, which the filter relies on):

```bash
seqtk seq -A -l 0 input.fasta | awk '/^>/{h=$0; next} {gsub(/[-.]/,"",$0); if (length($0)>0) print h"\n"$0}' > cleaned.fasta
```

## Getting Help

- [GitHub Issues](https://github.com/InstitutoTodosPelaSaude/viralQC/issues)
- Check logs in `output/logs/` directory
