# Developer Guide

This guide provides detailed documentation for the Python scripts used within the ViralQC pipeline. It is intended for developers who want to understand the internal logic, implementation choices, and execution flow of the tool.

The scripts are located in `viralqc/scripts/python/` and are orchestrated by Snakemake workflows.

## Dataset Management

These scripts are used by the `get-nextclade-datasets` and `get-blast-database` commands to download and prepare reference data.

### get_github_dataset.py

**Purpose**: Downloads specific dataset folders from a GitHub repository without using the GitHub API (to avoid rate limits).

**Execution**: Called by `get_public_datasets.smk` when processing `github` configured viruses.

**Implementation Details**:
- **No API Usage**: Instead of using the GitHub API, it downloads the repository archive via `https://codeload.github.com/.../zip/refs/heads/main`.
- **Selective Extraction**: It streams the zip file and only extracts files that match the requested `dataset-path`.
- **Structure Flattening**: It handles the stripping of the root folder (e.g., `repo-main/`) to place files directly in the target directory.

### jsonl_to_gff.py

**Purpose**: Converts NCBI Datasets JSONL annotation reports into GFF3 format, ensuring compatibility with Nextclade.

**Execution**: Called by `get_blast_database.smk` after downloading NCBI RefSeq data.

**Key Functions**:
- `clean_cds_name(cds_name)`: Sanitizes CDS names by removing special characters, truncating to 20 characters, and standardizing formatting. This is crucial because Nextclade can fail with complex or overly long gene names.
- `jsonl_to_gff(...)`:
  - **Validation**: Checks if CDS lengths are multiples of 3. If not, the accession is marked as invalid and excluded.
  - **Grouping**: Groups split CDS entries (e.g., joined genes) by name to create single gene features.
  - **Gene vs CDS**: If CDS data is missing, it attempts to create a gene feature using the full genome length (if divisible by 3).

### get_minimizer_index.py

**Purpose**: Generates a minimizer index JSON file from reference FASTA files, allowing Nextclade to map sequences to external (github-hosted) datasets.

**Execution**: Called by `get_public_datasets.smk` for GitHub-hosted datasets.

**Implementation Details**:
- **Origin**: This is a simplified adaptation of the `minimizer.py` script from the Nextclade project.
- **Customization**: The `fasta_read` function was modified to inject the dataset name into the sequence record annotations. This ensures that the generated index correctly associates sequence minimizers with their corresponding local dataset paths.

## Analysis Pipeline

These scripts are executed during the main `vqc run` workflow to process sequences, identify viruses, and assess quality.

### format_nextclade_sort.py

**Purpose**: Processes the output of `nextclade sort` to link identified datasets with their local file paths and identify sequences that didn't match any dataset.

**Execution**: Run immediately after `nextclade sort` in `run_analysis.smk`.

**Key Functions**:
- `map_datasets_to_local_paths(...)`: Reads the YAML configuration to build a mapping between remote dataset names (e.g., from Nextclade) and local storage paths.
- `format_nextclade_output(...)`: Merges "local" and "external" sort results. It adds a `localDataset` column pointing to the directory of the identified virus.
- `write_unmapped_sequences(...)`: Extracts sequences that have no assigned dataset (`localDataset` is NaN) and writes their names to `unmapped_sequences.txt` for subsequent BLAST analysis.

### blast_wrapper.py

**Purpose**: A wrapper around the `blastn` command to safely handle FASTA headers containing spaces.

**Execution**: Run by the `blast` rule in `run_analysis.smk` for sequences that were not identified by Nextclade.

**Implementation Details**:
- **Sanitization**: BLAST can truncate headers at the first space, leading to ID mismatches. This script checks for spaces in headers.
- **Renaming Flow**:
  1. If spaces are found, it generates a temporary FASTA where sequences are renamed to simple indices (1, 2, 3...).
  2. It saves a mapping file (`mapping.tsv`) linking indices to original headers.
  3. Runs BLAST with the renamed file.
  4. Restores the original headers in the BLAST output TSV using the mapping file.

### reorder_cds.py

**Purpose**: Reorders the `cdsCoverage` string in Nextclade's TSV output to match the gene order defined in the GFF file.

**Execution**: Run after every `nextclade run` execution.

**Logic**:
- Nextclade output genes in an alphabetical order and omit genes with zero coverage.
- This script reads the GFF to establish the canonical gene order (`start` position).
- It parses the existing `cdsCoverage` (format `Gene:Cov,...`), reorders it, and inserts `Gene:0.0` for any missing genes. This ensures consistent column ordering for downstream processing.

### post_process_nextclade.py

**Purpose**: The central aggregation script that combines results from Nextclade, BLAST, and generic analyses into a final report. It calculates categorical quality metrics (grades A-D) and produces the final output (TSV, CSV, or JSON).

**Execution**: The final step of the `post_process_nextclade` rule.

**Memory Management & Generators**:

This script is engineered to process massive datasets (e.g., millions of sequences) with a constant memory footprint for CSV/TSV outputs.

1.  **Lazy Loading with Generators**:
    The `format_dfs` function is implemented as a Python **Generator**. Instead of returning a list of all DataFrames (which would load all files into RAM), it `yields` one processed DataFrame at a time.
    
    *   **Logic**: It iterates through the list of input files. For each file, it reads the data, enriches it with metadata, optimizes types, yields it to the consumer, and **immediately** deletes the reference and forces garbage collection.

2.  **Streamed Writing**:
    The `write_combined_df` (and its helper `_write_csv_tsv_output`) consumes this generator. It iterates over the generator, writing each yielded chunk to disk immediately using `mode='a'` (append).
    
    *   **Result**: At any given moment, only the data from a single input file exists in memory.
    
    *   **JSON Limitation**: For JSON output (`_write_json_output`), the script *must* accumulate all data to form a valid JSON array. However, it still employs garbage collection to discard intermediate processing artifacts as soon as they are appended to the main list.

3.  **Explicit Garbage Collection**:
    Python's automatic GC might not trigger fast enough when dealing with large, tight loops of data loading. Explicit `del` combined with `gc.collect()` calls are strategically placed to ensure memory is released back to the OS before allocating the next chunk.

**Core Functions**:

*   `format_dfs(files, config_file, blast_metadata_df)`:
    The primary generator. It determines if a result file belongs to a known virus (with dataset configured in YAML) or is a generic run (Nextclade executed with references informed by BLAST analysis). It calls the appropriate processing logic (`_process_with_virus_info` or `_process_generic_run`) and yields the result.

*   `load_blast_metadata(metadata_path)`:
    Loads the BLAST database metadata and normalizes column names (e.g., `accession` -> `virus`) to ensure consistency with Nextclade outputs.

*   `optimize_dataframe_memory(df)`:
    Analyzes DataFrame columns. If a string column (like `virus`, `clade`, `dataset`) has a low cardinality (number of unique values < 50% of total rows), it converts the column to `category` dtype. This drastically reduces RAM usage.

*   `add_qualities(df, virus_info)`:
    Applies the quality scoring logic row-by-row. It invokes the helper `_compute_metrics_qualities` and then uses the `get_*_quality` functions to assign grades.

*   `add_coverages(df, virus_info)`:
    Parses and Formats the `cdsCoverage` string. It also calculates coverage for specific target regions and adds them as new columns (`targetRegionsCoverage`, `targetGeneCoverage`).

*   `format_sc2_clade(df, dataset_name)`:
    Contains specific logic for SARS-CoV-2. Since Nextclade outputs Pango lineages in a specific column (`Nextclade_pango`), this function maps it to the standard `clade` column for consistency.

*   `create_unmapped_df(unmapped_sequences, blast_results, blast_metadata_df)`:
    Handles sequences that failed both Nextclade identification and BLAST. It reads the raw `unmapped_sequences.txt` and creates a DataFrame labeled as "Unclassified".

*   `write_combined_df(df_iterator, output_file, output_format, ...)`:
    The main dispatcher. It takes the `format_dfs` generator (iterator) and directs it to the appropriate writing backend (CSV/TSV or JSON).

**Quality Control Functions**:

*   `get_genome_quality(scores)`:
    Aggregates individual metric scores into a final Genome Quality. It sums the scores (A=4, B=3, C=2, D=1) and normalizes to a 24-point scale to assign the final grade.
*   `get_target_regions_quality(...)`:
    Determines the quality of specific target regions. Logic: If the whole genome is A/B, return empty (implied good). Otherwise, calculate mean coverage of the target regions to assign a grade.
*   `get_cds_cov_quality(...)`:
    Checks every CDS against coverage thresholds to assign A/B/C/D grades per gene.
*   `get_missing_data_quality(coverage)`: Scored based on thresholds (0.9, 0.75, 0.5).
*   `get_private_mutations_quality(total, threshold)`: Scored based on deviation from a defined threshold.
*   `get_qc_quality(total)`: General scoring for count-metrics (0=A, 1=B, 2=C, >2=D).

### extract_target_regions.py

**Purpose**: Extracts the genomic coordinates of "good quality" regions for downstream usage (e.g., consensus generation or primer design).

**Execution**: Run by `extract_target_regions` rule after post-processing.

**Selection Logic**:

The `check_target_regions` function determines the best region to extract based on quality:

1. **Genome**: If the overall `genomeQuality` is A or B, the entire genome is selected.
2. **Target Region**: Else, if `targetRegionsQuality` is A or B, the specific target regions are selected.
3. **Target Gene**: Else, if `targetGeneQuality` is A or B, the target gene is selected.

**Coordinate Mapping**:
- Uses `get_regions` to look up the start/end coordinates of the selected feature (gene or full genome) in the GFF file.
- Outputs a BED file compatible with `seqtk subseq`.
