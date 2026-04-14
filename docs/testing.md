# Testing

`viralqc` uses [pytest](https://docs.pytest.org) for unit testing. All tests live in the `tests/` directory at the repository root.

## Setup

Install the package and dev extras in the `viralQC` micromamba environment:

```bash
micromamba run -n viralQC pip install -e ".[dev]"
```

## Running the tests

From the repository root:

```bash
# All tests (verbose, short tracebacks — configured in pyproject.toml)
micromamba run -n viralQC pytest

# Only one module
micromamba run -n viralQC pytest tests/test_ncbi_submission.py

# Run a specific test class or function
micromamba run -n viralQC pytest tests/test_ncbi_submission.py::TestFilterSequences
micromamba run -n viralQC pytest tests/test_ncbi_submission.py::TestFilterSequences::test_high_n_content_dropped
```

## Test structure

```
tests/
├── conftest.py              # Shared fixtures and file-builder helpers
├── test_ncbi_submission.py  # Unit tests for viralqc/core/ncbi_submission.py
├── test_run_analysis.py     # Unit tests for viralqc/core/run_analysis.py
└── test_scripts_python.py   # Unit tests for viralqc/scripts/python/
```

### `conftest.py`

Provides shared helpers and pytest fixtures:

| Symbol | Description |
|--------|-------------|
| `make_fasta(tmp_path, name, records)` | Write a minimal FASTA file |
| `make_results_tsv(tmp_path, rows)` | Write a minimal ViralQC results TSV |
| `make_metadata_csv(tmp_path, rows)` | Write a minimal metadata CSV |
| `tmp_fasta` (fixture) | Two-sequence FASTA (both pass quality filters) |
| `tmp_results` (fixture) | SARS-CoV-2 results TSV with two rows |
| `tmp_metadata` (fixture) | Standard-virus metadata CSV with two rows |

### `test_ncbi_submission.py`

Covers every public function in `viralqc/core/ncbi_submission.py`:

| Test class | Functions under test |
|---|---|
| `TestSanitizeSeqName` | `sanitize_seq_name` |
| `TestNFraction` | `_n_fraction` |
| `TestLoadResults` | `load_results` (TSV / CSV / JSON / invalid) |
| `TestLoadFasta` | `load_fasta` |
| `TestFilterSequences` | `filter_sequences` (N-content, length, custom thresholds) |
| `TestWriteFasta` | `write_fasta` (single file, batching, rename log) |
| `TestBatchPath` | `_batch_path` |
| `TestWriteHelpers` | `write_dropped_table`, `write_rename_log` |
| `TestCopyTbl` | `copy_tbl` |
| `TestIsPlainHeaderVirus` | `is_plain_header_virus` |
| `TestSerotypeFields` | `_serotype_fields` |
| `TestLoadMetadata` | `load_metadata` |
| `TestWriteSubmissionMetadata` | `write_submission_metadata` |

### `test_run_analysis.py`

Covers `viralqc/core/run_analysis.RunAnalysis`. Snakemake is **mocked** with `unittest.mock.patch` so no external tools are required:

| Test class | What is tested |
|---|---|
| `TestGetOutputFormat` | Valid / invalid output extensions |
| `TestRun` | Snakemake is called, config keys are correct, absolute paths, invalid format raises |

### `test_scripts_python.py`

Tests pure-Python helper functions from each script. I/O-heavy entry points (`main()`, `__main__`) are not tested here — they rely on integration tests with real data.

| Test class | Script / function |
|---|---|
| `TestValidateFasta` | `validate_fasta.validate_fasta_file` |
| `TestParseCdsCoverage` / `TestReorderCdsCoverage` / `TestReadGffGeneOrder` / `TestProcessNextcladeTsv` | `reorder_cds` |
| `TestSanitizeName` / `TestLoadIdMapping` / `TestParseTblBlocks` | `split_tbl_by_sample` |
| `TestCreateFastaPath` / `TestMapDatasetsToLocalPaths` / `TestWriteUnmappedSequences` | `format_nextclade_sort` |
| `TestParseFastaLengths` / `TestCleanCdsName` | `jsonl_to_gff` |

## Writing new tests

1. Place test files in `tests/` with the prefix `test_`.
2. Import helper utilities from `tests/conftest.py` when you need mock files.
3. Use `tmp_path` (built-in pytest fixture) for temporary directories — pytest cleans them up automatically.
4. Prefer parameterised tests (`@pytest.mark.parametrize`) for repetitive cases.
5. Mock external calls (`run_snakemake`, network requests, etc.) with `unittest.mock.patch`.

```python
# Example: mock an external dependency
from unittest.mock import patch, MagicMock

@patch("viralqc.core.run_analysis.run_snakemake")
def test_something(mock_snake, tmp_path):
    mock_snake.return_value = MagicMock(success=True)
    ...
```
