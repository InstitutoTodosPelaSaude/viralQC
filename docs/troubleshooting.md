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

**Symptom:** a transfer from NCBI is interrupted while building the BLAST database
```
Error: Download error: stream error: stream ID 3; INTERNAL_ERROR; received from peer
Error in rule makeblast_db:
ERROR: Failed to create BLAST database.
```

NCBI aborts transfers intermittently, which is more likely on a slow connection or
when the download runs on every deploy. `vqc get-blast-database` retries each NCBI
download up to 3 times, waiting 10 seconds between attempts, so a single interrupted
transfer no longer fails the run.

**Solutions:**
1. Re-run `vqc get-blast-database` — all 3 attempts failing usually means NCBI is
   temporarily unavailable rather than a problem with your setup
2. Check internet connection and any proxy or firewall between you and
   `ftp.ncbi.nih.gov`
3. If you rebuild containers frequently, keep the datasets directory on a persistent
   volume so it is downloaded once instead of on every start

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

**Cause:** `blastn`'s FASTA reader rejects the whole query file when a line looks
implausible, and it does so on a sliding scale: measured against a real database
it aborts once a line falls below roughly **one third** plausible residues — not
only when a record is 100% gaps. So both an empty consensus exported as an
all-gap alignment row and a heavily gapped partial genome could trigger it.
Because the rule runs under `set -euo pipefail`, that one record aborted the
whole run and discarded every other sequence's results.

**Solution:** upgrade to **v1.2.0 or later**, which addresses this in two layers:

1. Alignment gaps (`-` and `.`) are stripped from every sequence, and records
   left with no usable base are excluded from Nextclade and BLAST, listed in
   `excluded_sequences.tsv`, and reported in `results.tsv` with
   `inputSequenceStatus` set to `no usable sequence`.
2. A non-zero `blastn` exit is no longer fatal. It is logged to
   `output/logs/blast.log` and the run continues with no BLAST hits, so other
   padding conventions that ViralQC does not strip — NEXUS-style `?`, for
   instance — degrade to "unidentified" instead of discarding the whole run.
   Affected sequences appear as `Unclassified`.

No change to the input file is required.

On earlier versions, strip the gaps and drop the empty records before running
(`-l 0` puts each sequence on a single line, which the filter relies on):

```bash
seqtk seq -A -l 0 input.fasta | awk '/^>/{h=$0; next} {gsub(/[-.]/,"",$0); if (length($0)>0) print h"\n"$0}' > cleaned.fasta
```

## Getting Help

- [GitHub Issues](https://github.com/InstitutoTodosPelaSaude/viralQC/issues)
- Check logs in `output/logs/` directory
