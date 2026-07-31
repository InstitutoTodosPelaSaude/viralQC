import argparse, re, csv, os, gc
from pathlib import Path as _Path
from itertools import chain
from pathlib import Path
from typing import Generator, Iterator
from pandas import read_csv, concat, DataFrame, notna, Series, NA, to_numeric
from numpy import nan
from pandas.errors import EmptyDataError
from yaml import safe_load


def sanitize_name(name: str) -> str:
    """
    Replace characters that are unsafe in file names with underscores.

    Args:
        name: Original sample header string.

    Returns:
        String with special characters replaced by '_'.
    """
    return re.sub(r"[^A-Za-z0-9_.\-]", "_", name)


TARGET_COLUMNS = {
    "seqName": str,
    "virus": str,
    "virus_tax_id": "Int64",
    "virus_species": str,
    "virus_species_tax_id": "Int64",
    "segment": str,
    "ncbi_id": str,
    "clade": str,
    # Alternative clade/lineage designations emitted by various Nextclade datasets.
    # Retained so a build can select the right one via viralqc.clade_column; reindex
    # in _sanitize_dataframe fills NaN for datasets that do not emit a given column.
    # Motivation: Nextclade's 2025-10-22 seasonal-flu change made the new short
    # "subclade" the default `clade` and moved the classic 6B.1A.5a.2a to
    # `legacy-clade`. Column names vary by dataset: hMPV uses the underscore form
    # `legacy_clade`; the combined flu-B dataset uses per-lineage suffixes
    # (`legacy-clade-vic`/`-yam`); mpox/Ebola carry finer `lineage`/`outbreak`.
    "legacy-clade": str,
    "legacy_clade": str,
    "legacy-clade-vic": str,
    "legacy-clade-yam": str,
    "subclade": str,
    "short-clade": str,
    "proposedSubclade": str,
    "lineage": str,
    "outbreak": str,
    "clade_who": str,
    "clade_nextstrain": str,
    "clade_display": str,
    "targetRegions": str,
    "targetGene": str,
    "genomeQuality": str,
    "genomeQualityScore": str,
    "targetRegionsQuality": str,
    "targetGeneQuality": str,
    "cdsCoverageQuality": str,
    "missingDataQuality": str,
    "privateMutationsQuality": str,
    "mixedSitesQuality": str,
    "snpClustersQuality": str,
    "frameShiftsQuality": str,
    "stopCodonsQuality": str,
    "coverage": "float64",
    "cdsCoverage": str,
    "targetRegionsCoverage": str,
    "targetGeneCoverage": str,
    "qc.overallScore": "float64",
    "qc.overallStatus": str,
    "alignmentScore": "float64",
    "substitutions": str,
    "deletions": str,
    "insertions": str,
    "frameShifts": str,
    "aaSubstitutions": str,
    "aaDeletions": str,
    "aaInsertions": str,
    "totalSubstitutions": "Int64",
    "totalDeletions": "Int64",
    "totalInsertions": "Int64",
    "totalFrameShifts": "Int64",
    "totalMissing": "Int64",
    "totalNonACGTNs": "Int64",
    "totalAminoacidSubstitutions": "Int64",
    "totalAminoacidDeletions": "Int64",
    "totalAminoacidInsertions": "Int64",
    "totalUnknownAa": "Int64",
    "qc.privateMutations.total": "Int64",
    "privateNucMutations.totalLabeledSubstitutions": "Int64",
    "privateNucMutations.totalUnlabeledSubstitutions": "Int64",
    "privateNucMutations.totalReversionSubstitutions": "Int64",
    "privateNucMutations.totalPrivateSubstitutions": "Int64",
    "qc.privateMutations.score": "float64",
    "qc.privateMutations.status": str,
    "qc.missingData.score": "float64",
    "qc.missingData.status": str,
    "qc.mixedSites.totalMixedSites": "Int64",
    "qc.mixedSites.score": "float64",
    "qc.mixedSites.status": str,
    "qc.snpClusters.totalSNPs": "Int64",
    "qc.snpClusters.score": "float64",
    "qc.snpClusters.status": str,
    "qc.frameShifts.totalFrameShifts": "Int64",
    "qc.frameShifts.score": "float64",
    "qc.frameShifts.status": str,
    "qc.stopCodons.totalStopCodons": "Int64",
    "qc.stopCodons.score": "float64",
    "qc.stopCodons.status": str,
    "gff_path": str,
    "tbl_path": str,
    "dataset": str,
    "datasetVersion": str,
    # Appended last on purpose: consumers that read results.tsv by column
    # position keep working. Empty for every sequence that reached nextclade or
    # BLAST; set for sequences excluded upstream for holding no usable base.
    "inputSequenceStatus": str,
}

# Value written to inputSequenceStatus for sequences that itemize_sequences kept
# out of the analysis because they carry no usable base. Kept as a single stable
# string so results can be filtered programmatically; the per-sequence reason is
# in excluded_sequences.tsv.
UNUSABLE_SEQUENCE_STATUS = "no usable sequence"


DEFAULT_PRIVATE_MUTATION_TOTAL_THRESHOLD = 10
COVERAGES_THRESHOLD = {
    "A": 0.95,
    "B": 0.75,
    "C": 0.5,
}

# Columns that should use categorical type for memory efficiency
CATEGORICAL_COLUMNS = [
    "virus",
    "virus_species",
    "segment",
    "ncbi_id",
    "clade",
    "legacy-clade",
    "legacy_clade",
    "legacy-clade-vic",
    "legacy-clade-yam",
    "subclade",
    "short-clade",
    "proposedSubclade",
    "lineage",
    "outbreak",
    "clade_who",
    "clade_nextstrain",
    "clade_display",
    "dataset",
    "datasetVersion",
    "inputSequenceStatus",
]

# Clade/lineage designation columns. Read as strings so all-integer values
# (e.g. flu legacy clade 1, 6, 7) are never coerced to float ('1.0') at read
# time. Kept narrow (not dtype=str for the whole file) because numeric columns
# such as `coverage` are compared numerically before _sanitize_dataframe runs.
CLADE_LIKE_COLUMNS = (
    "clade",
    "legacy-clade",
    "legacy_clade",
    "legacy-clade-vic",
    "legacy-clade-yam",
    "subclade",
    "short-clade",
    "proposedSubclade",
    "lineage",
    "outbreak",
    "clade_who",
    "clade_nextstrain",
    "clade_display",
)


def load_blast_metadata(metadata_path: Path) -> DataFrame:
    """
    Load the BLAST metadata TSV file.

    Args:
        metadata_path: Path to the metadata file.

    Returns:
        Dataframe containing the BLAST metadata.
    """
    column_mapping = {
        "accession": "virus",
        "segment": "segment",
        "virus_name": "virus_name",
        "virus_tax_id": "virus_tax_id",
        "release_date": "release_date",
        "species_name": "species_name",
        "species_tax_id": "species_tax_id",
        "database_version": "dataset_with_version",
    }
    try:
        df = read_csv(metadata_path, sep="\t", header=0)
        df = df.rename(columns=column_mapping)
        df = df.loc[:, ~df.columns.duplicated()]

        return df
    except Exception:
        return DataFrame(columns=list(column_mapping.values()))


def load_id_mapping(mapping_path: Path) -> dict:
    """
    Load the ID mapping TSV file.

    Args:
        mapping_path: Path to the mapping file.

    Returns:
        Dictionary mapping sanitized ID to original header.
    """
    try:
        df = read_csv(mapping_path, sep="\t", dtype=str)
        return dict(zip(df["id"], df["original_header"]))
    except Exception:
        return {}


def format_sc2_clade(df: DataFrame, dataset_name: str) -> DataFrame:
    """
    For SARS-CoV-2 datasets, replace 'clade' with 'Nextclade_pango'.

    Args:
        df: Dataframe of nextclade results.
        dataset_name: Name of dataset.

    Returns:
        For SARS-CoV-2 datasets returns a dataframe with values from
        Nextclade_pango column into clade column.
    """
    if dataset_name.startswith("sarscov2"):
        if "Nextclade_pango" in df.columns:
            df["clade"] = df["Nextclade_pango"]

    return df


def get_missing_data_quality(coverage: float) -> str:
    """
    Calculate missing data quality score based on coverage.

    Args:
        coverage: Genome coverage value.

    Returns:
        Quality score ('A', 'B', 'C', 'D' or empty string).
    """
    if not notna(coverage):
        return ""
    elif coverage >= 0.9:
        return "A"
    elif coverage >= 0.75:
        return "B"
    elif coverage >= 0.5:
        return "C"
    else:
        return "D"


def get_private_mutations_quality(total: int, threshold: int) -> str:
    """
    Calculate private mutations quality score.

    Args:
        total: Total number of private mutations.
        threshold: Threshold for private mutations.

    Returns:
        Quality score ('A', 'B', 'C', 'D' or empty string).
    """
    if not notna(total):
        return ""
    elif total <= threshold:
        return "A"
    elif total <= threshold * 1.05:
        return "B"
    elif total <= threshold * 1.1:
        return "C"
    else:
        return "D"


def get_qc_quality(total: int) -> str:
    """
    Calculate general QC quality score based on total count.

    Args:
        total: Total count of the metric.

    Returns:
        Quality score ('A', 'B', 'C', 'D' or None).
    """
    if not notna(total):
        return None
    elif total == 0:
        return "A"
    elif total == 1:
        return "B"
    elif total == 2:
        return "C"
    else:
        return "D"


def get_genome_quality(scores: list[str]) -> tuple[int, str]:
    """
    Evaluate the quality of genome based on 6 quality scores.

    Args:
        scores: List of scores categories.

    Returns:
        The quality of the genome.
    """
    values = {"A": 4, "B": 3, "C": 2, "D": 1}
    valid_scores = [values[s] for s in scores if s in values]

    total = sum(valid_scores)

    if not valid_scores:
        return 0, ""

    max_possible = len(valid_scores) * 4
    normalized_total = (total / max_possible) * 24

    if normalized_total == 24:
        return normalized_total, "A"
    elif normalized_total >= 18:
        return normalized_total, "B"
    elif normalized_total >= 12:
        return normalized_total, "C"

    return normalized_total, "D"


def _parse_cds_cov(cds_list: str | dict) -> dict[str, float]:
    """
    Parse the cdsCoverage string into a dictionary.

    Args:
        cds_list: String or dict containing CDS coverage data.

    Returns:
        Dictionary mapping gene names to coverage values.
    """
    if isinstance(cds_list, dict):
        return cds_list
    if not isinstance(cds_list, str):
        return {}
    parts = cds_list.split(",")
    result = {}
    for p in parts:
        if ":" in p:
            cds, cov = p.split(":")
            try:
                result[cds.strip()] = round(float(cov), 4)
            except ValueError:
                continue
    return result


def _normalize_cds_coverage(cds_coverage: str | dict) -> dict[str, float]:
    """
    Normalize cds_coverage input to always return a dictionary.

    Args:
        cds_coverage: Value of the 'cdsCoverage' column (str or dict).

    Returns:
        Dictionary mapping gene names to coverage values.
    """
    if isinstance(cds_coverage, str):
        return _parse_cds_cov(cds_coverage)
    elif isinstance(cds_coverage, dict):
        return cds_coverage
    else:
        return {}


def get_cds_cov_quality(
    cds_coverage: str | dict,
    target_threshold_a: float,
    target_threshold_b: float,
    target_threshold_c: float,
) -> str:
    """
    Categorize the CDS regions based on coverage thresholds.

    Args:
        cds_coverage: Value of the 'cdsCoverage' column from the Nextclade output (str or dict).
        target_threshold_a: Minimum required coverage to consider a target region as "A".
        target_threshold_b: Minimum required coverage to consider a target region as "B".
        target_threshold_c: Minimum required coverage to consider a target region as "C".

    Returns:
        The status of the target regions.
    """
    cds_coverage = _normalize_cds_coverage(cds_coverage)
    if not cds_coverage:
        return ""

    result = {}
    for cds, cov in cds_coverage.items():
        try:
            cov_val = float(cov)
            if cov_val >= target_threshold_a:
                result[cds] = "A"
            elif cov_val >= target_threshold_b:
                result[cds] = "B"
            elif cov_val >= target_threshold_c:
                result[cds] = "C"
            elif cov_val > 0:
                result[cds] = "D"
        except (ValueError, TypeError):
            continue

    return ", ".join(f"{cds}: {coverage}" for cds, coverage in result.items())


def get_target_regions_quality(
    cds_coverage: str | dict,
    genome_quality: str,
    target_regions: list,
    target_threshold_a: float,
    target_threshold_b: float,
    target_threshold_c: float,
) -> str:
    """
    Evaluate the quality of target regions and classify them as categories based
    on coverage thresholds.

    Args:
        cds_coverage: Value of the 'cdsCoverage' column from the Nextclade output (str or dict).
        genome_quality: Quality of genome.
        target_regions: List of target regions.
        target_threshold_a: Minimum required coverage to consider a target region as "A".
        target_threshold_b: Minimum required coverage to consider a target region as "B".
        target_threshold_c: Minimum required coverage to consider a target region as "C".

    Returns:
        The status of the target regions.
    """
    if genome_quality in ["A", "B", ""]:
        return ""

    if not target_regions:
        return ""

    cds_coverage = _normalize_cds_coverage(cds_coverage)

    cds_coverage = {k.strip(): v for k, v in cds_coverage.items()}
    coverages = []
    for region in target_regions:
        coverages.append(float(cds_coverage.get(region, 0)))

    if not coverages:
        return ""

    mean_coverage = sum(coverages) / len(coverages)
    if mean_coverage >= target_threshold_a:
        return "A"
    elif mean_coverage >= target_threshold_b:
        return "B"
    elif mean_coverage >= target_threshold_c:
        return "C"

    return "D"


def get_target_regions_coverage(
    cds_coverage: str | dict, target_regions: list[str]
) -> str:
    """
    Extract the coverage of specific genomic regions.

    Args:
        cds_coverage: Value of the 'cdsCoverage' column from the Nextclade output (str or dict).
        target_regions: List of target regions.

    Returns:
        A string with region and coverage.
    """
    cds_coverage = _normalize_cds_coverage(cds_coverage)

    target_cds_coverage = [
        f"{region}: {cds_coverage.get(region,0)}" for region in target_regions
    ]

    return ", ".join(target_cds_coverage)


def add_coverages(df: DataFrame, virus_info: dict) -> None:
    """
    Add 'targetRegionsCoverage', 'targetGeneCoverage', and format
    'cdsCoverage' column to results dataframe (in-place).

    Args:
        df: Dataframe of nextclade results.
        virus_info: Dictionary with specific virus configuration
    """
    if "cdsCoverage" not in df.columns:
        df["cdsCoverage"] = ""

    df["targetRegionsCoverage"] = df["cdsCoverage"].apply(
        lambda cds_cov: (
            get_target_regions_coverage(cds_cov, virus_info["target_regions"])
            if notna(cds_cov)
            else ""
        )
    )

    target_gene = virus_info.get("target_gene")
    df["targetGeneCoverage"] = df["cdsCoverage"].apply(
        lambda cds_cov: (
            get_target_regions_coverage(cds_cov, [target_gene])
            if notna(cds_cov) and target_gene
            else ""
        )
    )

    # Format cdsCoverage as string (will be converted to array for JSON output later)
    df["cdsCoverage"] = df["cdsCoverage"].apply(_parse_cds_cov)
    df["cdsCoverage"] = df["cdsCoverage"].apply(
        lambda d: ", ".join(f"{cds}: {coverage}" for cds, coverage in d.items())
    )


def _compute_metrics_qualities(row: Series, virus_info: dict) -> dict[str, str]:
    """
    Compute individual QC metrics quality scores from a row.

    Args:
        row: DataFrame row containing QC metrics.
        virus_info: Virus configuration dictionary.

    Returns:
        Dictionary with quality scores for each metric.
    """
    missing_data_quality = get_missing_data_quality(row.get("coverage", nan))

    private_mutations_total = row.get("qc.privateMutations.total", nan)
    private_mutations_quality = get_private_mutations_quality(
        total=private_mutations_total,
        threshold=virus_info.get(
            "private_mutation_total_threshold",
            DEFAULT_PRIVATE_MUTATION_TOTAL_THRESHOLD,
        ),
    )

    mixed_sites_quality = get_qc_quality(row.get("qc.mixedSites.totalMixedSites", nan))
    snp_clusters_quality = get_qc_quality(row.get("qc.snpClusters.totalSNPs", nan))
    frameshifts_quality = get_qc_quality(
        row.get("qc.frameShifts.totalFrameShifts", nan)
    )
    stop_codons_quality = get_qc_quality(row.get("qc.stopCodons.totalStopCodons", nan))

    return {
        "missingDataQuality": missing_data_quality,
        "privateMutationsQuality": private_mutations_quality,
        "mixedSitesQuality": mixed_sites_quality,
        "snpClustersQuality": snp_clusters_quality,
        "frameShiftsQuality": frameshifts_quality,
        "stopCodonsQuality": stop_codons_quality,
    }


def _compute_genome_quality(metrics_qualities: dict) -> tuple[str, str]:
    """
    Compute overall genome quality from individual metrics.

    Args:
        metrics_qualities: Dictionary of individual quality scores.

    Returns:
        Tuple of (genome_score, genome_quality).
    """
    genome_score, genome_quality = get_genome_quality(
        [
            metrics_qualities["missingDataQuality"],
            metrics_qualities["mixedSitesQuality"],
            metrics_qualities["privateMutationsQuality"],
            metrics_qualities["snpClustersQuality"],
            metrics_qualities["frameShiftsQuality"],
            metrics_qualities["stopCodonsQuality"],
        ]
    )
    return genome_score, genome_quality


def _compute_target_qualities(
    row: Series, virus_info: dict, genome_quality: str
) -> dict[str, str]:
    """
    Compute target region and gene quality scores.

    Args:
        row: DataFrame row containing CDS coverage data.
        virus_info: Virus configuration dictionary.
        genome_quality: Overall genome quality score.

    Returns:
        Dictionary with target region quality scores.
    """
    cds_coverage = row.get("cdsCoverage", nan)

    if not notna(cds_coverage) or cds_coverage == "":
        return {
            "targetRegionsQuality": "",
            "targetGeneQuality": "",
            "cdsCoverageQuality": "",
        }

    target_regions_quality = get_target_regions_quality(
        cds_coverage=cds_coverage,
        genome_quality=genome_quality,
        target_regions=virus_info["target_regions"],
        target_threshold_a=COVERAGES_THRESHOLD["A"],
        target_threshold_b=COVERAGES_THRESHOLD["B"],
        target_threshold_c=COVERAGES_THRESHOLD["C"],
    )

    target_gene_quality = get_target_regions_quality(
        cds_coverage=cds_coverage,
        genome_quality=target_regions_quality,
        target_regions=[virus_info["target_gene"]],
        target_threshold_a=COVERAGES_THRESHOLD["A"],
        target_threshold_b=COVERAGES_THRESHOLD["B"],
        target_threshold_c=COVERAGES_THRESHOLD["C"],
    )

    cds_cov_quality = get_cds_cov_quality(
        cds_coverage=cds_coverage,
        target_threshold_a=virus_info.get("target_regions_cov", COVERAGES_THRESHOLD)[
            "A"
        ],
        target_threshold_b=virus_info.get("target_regions_cov", COVERAGES_THRESHOLD)[
            "B"
        ],
        target_threshold_c=virus_info.get("target_regions_cov", COVERAGES_THRESHOLD)[
            "C"
        ],
    )

    return {
        "targetRegionsQuality": target_regions_quality,
        "targetGeneQuality": target_gene_quality,
        "cdsCoverageQuality": cds_cov_quality,
    }


def add_qualities(df: DataFrame, virus_info: dict) -> None:
    """
    Compute all quality metrics and add to dataframe (in-place).

    Args:
        df: Dataframe of nextclade results.
        virus_info: Dictionary with specific virus configuration
    """

    def compute_all_qualities(row):
        metrics = _compute_metrics_qualities(row, virus_info)
        genome_score, genome_quality = _compute_genome_quality(metrics)
        targets = _compute_target_qualities(row, virus_info, genome_quality)

        return Series(
            {
                **metrics,
                "genomeQualityScore": genome_score,
                "genomeQuality": genome_quality,
                **targets,
            }
        )

    qualities_df = df.apply(compute_all_qualities, axis=1)

    for col in qualities_df.columns:
        df[col] = qualities_df[col]


def optimize_dataframe_memory(df: DataFrame) -> None:
    """
    Optimize DataFrame memory usage by converting appropriate columns to categorical.
    Modifies DataFrame in-place.

    Args:
        df: DataFrame to optimize
    """
    for col in CATEGORICAL_COLUMNS:
        if col in df.columns and df[col].dtype == "object":
            # Only convert to categorical if there are repeated values
            n_unique = df[col].nunique()
            n_total = len(df)
            if n_unique < n_total * 0.5:
                df[col] = df[col].astype("category")


def _extract_dataset_name(file: str) -> str:
    """
    Extract the dataset name from a nextclade result filename.

    Args:
        file: Path to nextclade output file.

    Returns:
        Dataset name (e.g., 'flu-h1n1-ha' or 'NC_000001.1').
    """
    dataset = re.sub("\.nextclade.tsv", "", re.sub(".*/", "", file))
    dataset = re.sub("\.generic", "", dataset)
    return dataset


def _process_with_virus_info(
    df: DataFrame, virus_dataset: str, virus_info: dict
) -> None:
    """
    Process dataframe with known virus configuration (in-place).

    Args:
        df: Dataframe to process.
        virus_dataset: Name of the virus dataset.
        virus_info: Configuration dictionary for this virus.
    """
    df = format_sc2_clade(df, virus_dataset)
    df["virus"] = virus_info["virus_name"]
    df["virus_tax_id"] = virus_info["virus_tax_id"]
    df["virus_species"] = virus_info["virus_species"]
    df["virus_species_tax_id"] = virus_info["virus_species_tax_id"]
    df["segment"] = virus_info["segment"]
    df["ncbi_id"] = virus_info["ncbi_id"]
    df["dataset"] = virus_info["dataset"]
    df["datasetVersion"] = virus_info["tag"]
    df["targetGene"] = virus_info["target_gene"]
    df["targetRegions"] = "|".join(virus_info["target_regions"])
    add_coverages(df, virus_info)
    add_qualities(df, virus_info)


def _fill_missing_columns(df: DataFrame, column_types: dict) -> None:
    """
    Add missing columns to dataframe with appropriate default values (in-place).

    Args:
        df: Dataframe to fill.
        column_types: Dictionary mapping column names to their types.
    """
    for col, dtype in column_types.items():
        if col not in df.columns:
            if dtype == str:
                df[col] = ""
            elif dtype in ("float64", "Int64", bool):
                df[col] = None
            else:
                df[col] = ""


def _process_generic_run(
    df: DataFrame, virus_dataset: str, blast_metadata_df: DataFrame = None
) -> None:
    """
    Process dataframe from generic nextclade run (in-place).

    Args:
        df: Dataframe to process.
        virus_dataset: Accession ID (e.g., 'NC_000001.1').
        blast_metadata_df: Optional BLAST metadata for enrichment.
    """
    df["ncbi_id"] = virus_dataset

    # Enrich with BLAST metadata if available
    if blast_metadata_df is not None:
        meta = blast_metadata_df[blast_metadata_df["virus"] == virus_dataset]
        if not meta.empty:
            row = meta.iloc[0]
            df["virus"] = row["virus_name"]
            df["virus_tax_id"] = row["virus_tax_id"]
            df["virus_species"] = row["species_name"]
            df["virus_species_tax_id"] = row["species_tax_id"]
            df["segment"] = row["segment"]
            df["dataset"] = row["dataset_with_version"].split("_")[0]
            df["datasetVersion"] = row["dataset_with_version"].split("_")[1]
        else:
            df["virus"] = virus_dataset
    else:
        df["virus"] = virus_dataset

    # Ensure segment is set
    if "segment" not in df.columns:
        df["segment"] = "Unsegmented"
    else:
        df["segment"] = df["segment"].fillna("Unsegmented")

    # Set generic run columns to nan/empty
    df["clade"] = nan
    df["qc.overallScore"] = nan
    df["qc.overallStatus"] = nan

    # Fill missing columns
    _fill_missing_columns(df, TARGET_COLUMNS)

    # Format cdsCoverage without calculating qualities
    mock_virus_info = {"target_regions": [], "target_gene": ""}
    add_coverages(df, mock_virus_info)

    # Set quality columns to empty/None
    quality_cols = [
        "missingDataQuality",
        "privateMutationsQuality",
        "mixedSitesQuality",
        "snpClustersQuality",
        "frameShiftsQuality",
        "stopCodonsQuality",
        "genomeQualityScore",
        "genomeQuality",
        "targetRegionsQuality",
        "targetGeneQuality",
        "cdsCoverageQuality",
    ]
    for col in quality_cols:
        if col in TARGET_COLUMNS:
            if TARGET_COLUMNS[col] == str:
                df[col] = ""
            else:
                df[col] = None


def format_dfs(
    files: list[str], config_file: Path, blast_metadata_df: DataFrame = None
) -> Generator[DataFrame, None, None]:
    """
    Load and format nextclade outputs based on information defined for each virus.

    Args:
        files: List of paths of nextclade outputs.
        config_file: Path to the YAML configuration file listing nextclade datasets.
        blast_metadata_df: Dataframe with BLAST metadata (optional).

    Yields:
        Formatted dataframes one at a time for memory efficiency.
    """
    with config_file.open("r") as f:
        config = safe_load(f)

    for file in files:
        try:
            # Force clade-like columns to string so an all-integer legacy-clade
            # (flu `1`,`6`,`7`) is not inferred as float and rendered `1.0`,`6.0`
            # (which would split into phantom clade categories). Other columns keep
            # inferred dtypes — `coverage` etc. are used numerically before sanitize.
            df = read_csv(
                file, sep="\t", header=0, dtype={c: str for c in CLADE_LIKE_COLUMNS}
            )
            df = df.loc[:, ~df.columns.duplicated()]
        except EmptyDataError:
            df = DataFrame(columns=list(TARGET_COLUMNS.keys()))

        if not df.empty:
            virus_dataset = _extract_dataset_name(file)
            virus_info = config["nextclade_data"].get(virus_dataset)
            if virus_info is None and "github" in config:
                virus_info = config["github"].get(virus_dataset)

            if virus_info:
                _process_with_virus_info(df, virus_dataset, virus_info)
            else:
                _process_generic_run(df, virus_dataset, blast_metadata_df)

        df = df.loc[:, ~df.columns.duplicated()]
        optimize_dataframe_memory(df)

        yield df

        del df
        gc.collect()


def create_unmapped_df(
    unmapped_sequences: Path, blast_results: Path, blast_metadata_df: DataFrame
) -> DataFrame:
    """
    Create a dataframe of unmapped sequences.

    Args:
        unmapped_sequences: Path to unmapped_sequences.txt file.
        blast_results: Path to blast results of unmapped_sequences.txt.
        blast_metadata_df: Dataframe with BLAST metadata.

    Returns:
        A dataframe of unmapped sequences.
    """
    data = []
    with open(unmapped_sequences, "r") as f:
        for line in f:
            data.append((line.strip().strip('"').strip("'"), "Unclassified"))

    df = DataFrame(data, columns=["seqName", "virus"])
    del data
    gc.collect()

    for col in TARGET_COLUMNS.keys():
        if col not in df.columns:
            if TARGET_COLUMNS[col] == str:
                df[col] = ""
            elif TARGET_COLUMNS[col] == "float64":
                df[col] = None
            elif TARGET_COLUMNS[col] == "Int64":
                df[col] = None
            elif TARGET_COLUMNS[col] == bool:
                df[col] = None
            else:
                df[col] = ""

    if os.path.getsize(blast_results) == 0:
        optimize_dataframe_memory(df)
        return df.loc[:, ~df.columns.duplicated()]
    else:
        blast_columns = [
            "seqName",
            "qlen",
            "virus",
            "slen",
            "qstart",
            "qend",
            "sstart",
            "send",
            "evalue",
            "bitscore",
            "pident",
            "qcovs",
            "qcovhsp",
        ]

        blast_df = read_csv(blast_results, sep="\t", header=None, names=blast_columns)
        blast_df = blast_df.loc[:, ~blast_df.columns.duplicated()]

        # Use the passed metadata dataframe
        blast_df = blast_df.merge(blast_metadata_df, on="virus", how="left")
        blast_df = blast_df[
            [
                "seqName",
                "virus",
                "segment",
                "virus_name",
                "virus_tax_id",
                "species_name",
                "species_tax_id",
                "dataset_with_version",
            ]
        ]

        df["seqName"] = df["seqName"].astype(str)
        blast_df["seqName"] = blast_df["seqName"].astype(str)
        merged = df.merge(blast_df, on="seqName", how="left", suffixes=("_df1", "_df2"))

        del df, blast_df
        gc.collect()

        merged["virus"] = merged["virus_df2"].fillna(merged["virus_df1"])

        final_df = merged.drop(columns=["virus_df1", "virus_df2"])
        del merged
        gc.collect()

        final_df = final_df.assign(
            ncbi_id=final_df["virus"],
            virus=final_df["virus_name"].fillna("Unclassified").astype(str),
            virus_tax_id=final_df["virus_tax_id_df2"].astype("Int64"),
            virus_species=final_df["species_name"].fillna("Unclassified").astype(str),
            virus_species_tax_id=final_df["species_tax_id"].astype("Int64"),
            segment=final_df["segment_df2"].fillna("Unsegmented").astype(str),
        )
        split_result = final_df["dataset_with_version"].str.split("_", n=1, expand=True)
        if split_result.shape[1] == 2:
            final_df[["dataset", "datasetVersion"]] = split_result
        else:
            final_df["dataset"] = None
            final_df["datasetVersion"] = None

        final_df = final_df.loc[:, ~final_df.columns.duplicated()]
        optimize_dataframe_memory(final_df)

        return final_df


def create_excluded_df(excluded_sequences: Path) -> DataFrame:
    """
    Create a dataframe of sequences excluded from the analysis for holding no
    usable base.

    These sequences never reach nextclade sort, so they appear in neither the
    nextclade results nor unmapped_sequences.txt. Without this dataframe they
    would be missing from the results file entirely, silently reducing the row
    count below the number of input sequences.

    Args:
        excluded_sequences: Path to the excluded_sequences.tsv file written by
            itemize_sequences.

    Returns:
        A dataframe with one row per excluded sequence, flagged as a
        quality-control failure. Empty when there is nothing to report.
    """
    # No path given: the caller is not tracking excluded sequences at all.
    if not excluded_sequences:
        return DataFrame()

    # A path WAS given, so anything unreadable or malformed is abnormal. Failing
    # loudly matches the standard _sanitize_dataframe sets for a missing required
    # column: silently returning an empty frame here would drop every excluded
    # sequence from the results with no trace, which is the exact failure this
    # function exists to prevent.
    if not os.path.exists(excluded_sequences):
        raise FileNotFoundError(
            f"ViralQC post-processing: excluded sequences file {excluded_sequences} "
            "was passed but does not exist; sequences excluded for holding no "
            "usable base would be missing from the results."
        )

    try:
        excluded = read_csv(excluded_sequences, sep="\t", dtype=str)
    except EmptyDataError as error:
        raise ValueError(
            f"ViralQC post-processing: excluded sequences file {excluded_sequences} "
            "is empty and carries no header; itemize_sequences always writes one, "
            "so the file is truncated or corrupted."
        ) from error

    if "id" not in excluded.columns:
        raise KeyError(
            f"ViralQC post-processing: excluded sequences file {excluded_sequences} "
            f"has no 'id' column (found {list(excluded.columns)}); it is corrupted."
        )

    # Header present with no data rows is the normal case: nothing was excluded.
    if excluded.empty:
        return DataFrame()

    # seqName carries the sanitized numeric id so the frame flows through the
    # same _map_ids step as every other frame and is written under the original
    # header.
    df = DataFrame({"seqName": excluded["id"].astype(str)})

    _fill_missing_columns(df, TARGET_COLUMNS)

    df["virus"] = "Unclassified"
    df["virus_species"] = "Unclassified"
    df["segment"] = "Unsegmented"
    df["coverage"] = 0.0
    df["genomeQuality"] = "D"
    df["inputSequenceStatus"] = UNUSABLE_SEQUENCE_STATUS

    optimize_dataframe_memory(df)
    return df.loc[:, ~df.columns.duplicated()]


def _sanitize_dataframe(df: DataFrame) -> DataFrame:
    """
    Sanitize and prepare dataframe for output (removes duplicates, casts types, cleans data).

    Args:
        df: Dataframe to sanitize.

    Returns:
        Sanitized dataframe.
    """
    df = df.loc[:, ~df.columns.duplicated()]

    for col, dtype in TARGET_COLUMNS.items():
        if col in df.columns:
            if dtype in ("Int64", "float64"):
                df[col] = to_numeric(df[col], errors="coerce").astype(dtype)
            elif dtype == str:
                df[col] = df[col].astype("string")

    # reindex (not df[[...]]) so datasets whose Nextclade output lacks an
    # optional column (e.g. 'legacy-clade', present only for influenza HA)
    # get it filled with NaN instead of raising KeyError.
    #
    # reindex silently tolerates ANY missing column, which would mask a genuine
    # upstream schema break. Guard the small set of columns that every processed
    # chunk must carry so such a break fails loudly instead of emitting all-null.
    _REQUIRED = ("seqName", "clade", "coverage", "genomeQuality")
    _missing_required = [c for c in _REQUIRED if c not in df.columns]
    if _missing_required:
        raise KeyError(
            "ViralQC post-processing: required column(s) "
            f"{_missing_required} absent from the processed Nextclade output; "
            "the dataset schema may have changed."
        )
    final_df = df.reindex(columns=list(TARGET_COLUMNS.keys())).round(4)
    final_df = final_df.replace(r"^\s*$", nan, regex=True)
    return final_df


def _map_ids(df: DataFrame, id_map: dict) -> DataFrame:
    """
    Restore original IDs using the mapping dictionary.
    Modifies DataFrame in-place.

    Args:
        df: DataFrame to process.
        id_map: Dictionary mapping sanitized ID to original header.

    Returns:
        DataFrame with restored IDs.
    """
    if id_map and "seqName" in df.columns:
        # Convert seqName to string to ensure matching works
        df["seqName"] = df["seqName"].astype(str)
        df["seqName"] = df["seqName"].map(id_map).fillna(df["seqName"])
    return df


def _filter_unmapped_sequences(
    unmapped_df: DataFrame, processed_seq_names: set
) -> DataFrame:
    """
    Filter unmapped sequences to exclude already processed sequences.

    Args:
        unmapped_df: Unmapped sequences dataframe.
        processed_seq_names: Set of sequence names already processed.

    Returns:
        Filtered unmapped dataframe.
    """
    if processed_seq_names and "seqName" in unmapped_df.columns:
        return unmapped_df[
            ~unmapped_df["seqName"].astype(str).isin(processed_seq_names)
        ]
    return unmapped_df


def add_annotation_paths(
    df: DataFrame,
    reverse_id_map: dict,
    gff_per_sample_dir: Path | None,
    tbl_per_sample_dir: Path | None,
) -> None:
    """
    Add gff_path and tbl_path columns to a dataframe (in-place).

    Paths are constructed as ``<dir>/{id}_{sanitized_seqName}.{ext}`` and
    set to an empty string when the corresponding file does not exist or
    when no annotation directory was provided.

    Args:
        df: DataFrame whose seqName column already holds original sample
            headers (i.e. after _map_ids has been called).
        reverse_id_map: Dictionary mapping original_header to numeric id.
        gff_per_sample_dir: Path to the gff_files/per_sample/ directory.
        tbl_per_sample_dir: Path to the tbl_files/per_sample/ directory.
    """

    def _build_path(seq_name: str, directory: Path | None, ext: str) -> str:
        if directory is None or not isinstance(seq_name, str):
            return ""
        sample_id = reverse_id_map.get(seq_name, "")
        if not sample_id:
            return ""
        candidate = directory / f"{sample_id}_{sanitize_name(seq_name)}.{ext}"
        return str(candidate) if candidate.exists() else ""

    if "seqName" not in df.columns:
        df["gff_path"] = ""
        df["tbl_path"] = ""
        return

    df["gff_path"] = (
        df["seqName"]
        .astype(str)
        .apply(lambda s: _build_path(s, gff_per_sample_dir, "gff"))
    )
    df["tbl_path"] = (
        df["seqName"]
        .astype(str)
        .apply(lambda s: _build_path(s, tbl_per_sample_dir, "tbl"))
    )


def _format_json_columns(df: DataFrame) -> None:
    """
    Format specific columns for JSON output (convert strings to arrays/dicts).
    Modifies dataframe in-place.

    Args:
        df: Dataframe to format.
    """
    coverage_cols = ["cdsCoverageQuality", "cdsCoverage", "targetRegionsCoverage"]
    mutation_cols = [
        "substitutions",
        "deletions",
        "insertions",
        "frameShifts",
        "aaSubstitutions",
        "aaDeletions",
        "aaInsertions",
    ]

    for col in coverage_cols:
        if col in df.columns:
            if col == "cdsCoverage":
                df[col] = df[col].apply(
                    lambda val: (
                        [{k: v} for k, v in _parse_cds_cov(val).items()]
                        if isinstance(val, str) and val.strip()
                        else None
                    )
                )
            else:

                def parse_coverage_to_dicts(val):
                    if not isinstance(val, str) or not val.strip():
                        return None
                    result = []
                    for item in val.split(","):
                        item = item.strip()
                        if ":" in item:
                            parts = item.split(":", 1)
                            region = parts[0].strip()
                            try:
                                value = float(parts[1].strip())
                                result.append({region: value})
                            except ValueError:
                                value = parts[1].strip()
                                result.append({region: value})
                    return result if result else None

                df[col] = df[col].apply(parse_coverage_to_dicts)

    for col in mutation_cols:
        if col in df.columns:
            df[col] = df[col].apply(
                lambda val: (
                    val.split(",") if isinstance(val, str) and val.strip() else None
                )
            )


def _finalize_synthesised_frame(
    df: DataFrame,
    id_map: dict = None,
    reverse_id_map: dict = None,
    gff_per_sample_dir: Path | None = None,
    tbl_per_sample_dir: Path | None = None,
) -> DataFrame:
    """
    Apply the ID restoration, annotation-path and schema handling used for
    nextclade frames to a frame synthesised by ViralQC (unmapped or excluded
    sequences).

    Args:
        df: Frame whose seqName column still holds sanitized numeric ids.
        id_map: Dictionary for ID mapping (numeric id -> original header).
        reverse_id_map: Dictionary mapping original header to numeric id.
        gff_per_sample_dir: Path to the gff_files/per_sample/ directory.
        tbl_per_sample_dir: Path to the tbl_files/per_sample/ directory.

    Returns:
        The sanitized frame, ready to be written.
    """
    if id_map:
        df = _map_ids(df, id_map)
    add_annotation_paths(
        df, reverse_id_map or {}, gff_per_sample_dir, tbl_per_sample_dir
    )
    return _sanitize_dataframe(df)


def _empty_results_frame() -> DataFrame:
    """
    Build an empty frame carrying the full output schema.

    Used so a run whose sequences all failed still produces a readable results
    file with a header rather than no file at all.

    Returns:
        An empty dataframe with every target column.
    """
    return DataFrame(columns=list(TARGET_COLUMNS.keys()))


def _write_json_output(
    df_iterator: Iterator[DataFrame],
    unmapped_df: DataFrame,
    processed_seq_names: set,
    output_file: Path,
    id_map: dict = None,
    reverse_id_map: dict = None,
    gff_per_sample_dir: Path | None = None,
    tbl_per_sample_dir: Path | None = None,
    excluded_df: DataFrame = None,
) -> None:
    """
    Write all dataframes to JSON output file.

    Args:
        df_iterator: Iterator of dataframes to write.
        unmapped_df: Optional unmapped sequences dataframe.
        processed_seq_names: Set of already processed sequence names.
        output_file: Path to output file.
        id_map: Dictionary for ID mapping (numeric id -> original header).
        reverse_id_map: Dictionary mapping original header to numeric id.
        gff_per_sample_dir: Path to the gff_files/per_sample/ directory.
        tbl_per_sample_dir: Path to the tbl_files/per_sample/ directory.
        excluded_df: Optional dataframe of sequences excluded for holding no
            usable base.
    """
    all_data = []

    for df in df_iterator:
        if id_map:
            df = _map_ids(df, id_map)
        add_annotation_paths(
            df, reverse_id_map or {}, gff_per_sample_dir, tbl_per_sample_dir
        )
        final_df = _sanitize_dataframe(df)
        all_data.append(final_df)
        del df, final_df
        gc.collect()

    if unmapped_df is not None and not unmapped_df.empty:
        unmapped_df = _filter_unmapped_sequences(unmapped_df, processed_seq_names)
        if not unmapped_df.empty:
            all_data.append(
                _finalize_synthesised_frame(
                    unmapped_df,
                    id_map,
                    reverse_id_map,
                    gff_per_sample_dir,
                    tbl_per_sample_dir,
                )
            )

    if excluded_df is not None and not excluded_df.empty:
        # Cannot collide inside the workflow (excluded records never enter the
        # sanitized FASTA, so nextclade never sees them), but this script is also a
        # standalone CLI where a stale --excluded-sequences file would duplicate rows.
        excluded_df = _filter_unmapped_sequences(excluded_df, processed_seq_names)
    if excluded_df is not None and not excluded_df.empty:
        all_data.append(
            _finalize_synthesised_frame(
                excluded_df,
                id_map,
                reverse_id_map,
                gff_per_sample_dir,
                tbl_per_sample_dir,
            )
        )

    if not all_data:
        all_data.append(_empty_results_frame())

    combined_df = concat(all_data, ignore_index=True)
    combined_df = combined_df.loc[:, ~combined_df.columns.duplicated()]
    combined_df = combined_df.sort_values(by=["virus"])

    del all_data
    gc.collect()

    _format_json_columns(combined_df)

    json_content = combined_df.to_json(orient="table", indent=4)
    json_content = json_content.replace("\\/", "/")
    with open(output_file, "w", encoding="utf-8") as f:
        f.write(json_content)

    del combined_df
    gc.collect()


def _write_csv_tsv_chunk(
    df: DataFrame,
    output_file: Path,
    output_format: str,
    is_first: bool,
    is_header: bool,
) -> None:
    """
    Write a single dataframe chunk to CSV/TSV file.

    Args:
        df: Dataframe to write.
        output_file: Path to output file.
        output_format: 'csv' or 'tsv'.
        is_first: Whether this is the first chunk (affects write mode).
        is_header: Whether to write header row.
    """
    if output_format == "tsv":
        df.to_csv(
            output_file,
            sep="\t",
            index=False,
            header=is_header,
            mode="w" if is_first else "a",
        )
    elif output_format == "csv":
        df.to_csv(
            output_file,
            sep=";",
            index=False,
            header=is_header,
            mode="w" if is_first else "a",
            quoting=csv.QUOTE_NONNUMERIC,
        )


def _write_csv_tsv_output(
    df_iterator: Iterator[DataFrame],
    unmapped_df: DataFrame,
    processed_seq_names: set,
    output_file: Path,
    output_format: str,
    id_map: dict = None,
    reverse_id_map: dict = None,
    gff_per_sample_dir: Path | None = None,
    tbl_per_sample_dir: Path | None = None,
    excluded_df: DataFrame = None,
) -> None:
    """
    Write all dataframes to CSV/TSV output file incrementally.

    Args:
        df_iterator: Iterator of dataframes to write.
        unmapped_df: Optional unmapped sequences dataframe.
        processed_seq_names: Set of already processed sequence names.
        output_file: Path to output file.
        output_format: 'csv' or 'tsv'.
        id_map: Dictionary for ID mapping (numeric id -> original header).
        reverse_id_map: Dictionary mapping original header to numeric id.
        gff_per_sample_dir: Path to the gff_files/per_sample/ directory.
        tbl_per_sample_dir: Path to the tbl_files/per_sample/ directory.
        excluded_df: Optional dataframe of sequences excluded for holding no
            usable base.
    """
    first_chunk = True

    for df in df_iterator:
        if id_map:
            df = _map_ids(df, id_map)
        add_annotation_paths(
            df, reverse_id_map or {}, gff_per_sample_dir, tbl_per_sample_dir
        )
        final_df = _sanitize_dataframe(df)
        _write_csv_tsv_chunk(
            final_df, output_file, output_format, first_chunk, first_chunk
        )
        first_chunk = False
        del df, final_df
        gc.collect()

    if unmapped_df is not None and not unmapped_df.empty:
        unmapped_df = _filter_unmapped_sequences(unmapped_df, processed_seq_names)
        if not unmapped_df.empty:
            final_unmapped = _finalize_synthesised_frame(
                unmapped_df,
                id_map,
                reverse_id_map,
                gff_per_sample_dir,
                tbl_per_sample_dir,
            )
            # first_chunk, not a hardcoded False: when no nextclade frame was
            # written this block owns the header and must create the file.
            _write_csv_tsv_chunk(
                final_unmapped, output_file, output_format, first_chunk, first_chunk
            )
            first_chunk = False
            del final_unmapped
            gc.collect()

    if excluded_df is not None and not excluded_df.empty:
        # See the note in _write_json_output: guards the standalone-CLI case.
        excluded_df = _filter_unmapped_sequences(excluded_df, processed_seq_names)
    if excluded_df is not None and not excluded_df.empty:
        final_excluded = _finalize_synthesised_frame(
            excluded_df,
            id_map,
            reverse_id_map,
            gff_per_sample_dir,
            tbl_per_sample_dir,
        )
        _write_csv_tsv_chunk(
            final_excluded, output_file, output_format, first_chunk, first_chunk
        )
        first_chunk = False
        del final_excluded
        gc.collect()

    if first_chunk:
        # Nothing at all to report: still emit a header-only file so the run
        # ends with a readable result instead of a missing-output failure.
        _write_csv_tsv_chunk(
            _empty_results_frame(), output_file, output_format, True, True
        )


def write_combined_df(
    df_iterator: Iterator[DataFrame],
    output_file: Path,
    output_format: str,
    unmapped_df: DataFrame = None,
    processed_seq_names: set = None,
    id_map: dict = None,
    reverse_id_map: dict = None,
    gff_per_sample_dir: Path | None = None,
    tbl_per_sample_dir: Path | None = None,
    excluded_df: DataFrame = None,
) -> None:
    """
    Write dataframes incrementally to reduce memory usage.

    Args:
        df_iterator: Iterator of formatted dataframes.
        output_file: Path to output file.
        output_format: Format to write output (csv, tsv, or json).
        unmapped_df: Optional unmapped sequences dataframe.
        processed_seq_names: Set of sequence names already processed.
        id_map: Dictionary for ID mapping (numeric id -> original header).
        reverse_id_map: Dictionary mapping original header to numeric id.
        gff_per_sample_dir: Path to the gff_files/per_sample/ directory.
        tbl_per_sample_dir: Path to the tbl_files/per_sample/ directory.
        excluded_df: Optional dataframe of sequences excluded for holding no
            usable base.
    """
    if output_format == "json":
        _write_json_output(
            df_iterator,
            unmapped_df,
            processed_seq_names,
            output_file,
            id_map,
            reverse_id_map,
            gff_per_sample_dir,
            tbl_per_sample_dir,
            excluded_df=excluded_df,
        )
    else:
        _write_csv_tsv_output(
            df_iterator,
            unmapped_df,
            processed_seq_names,
            output_file,
            output_format,
            id_map,
            reverse_id_map,
            gff_per_sample_dir,
            tbl_per_sample_dir,
            excluded_df=excluded_df,
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process Nextclade output files.")

    parser.add_argument(
        "--files", nargs="*", default=[], help="List of Nextclade output .tsv files"
    )
    parser.add_argument(
        "--generic-files",
        nargs="*",
        default=[],
        help="List of Generic Nextclade output .tsv files",
    )
    parser.add_argument(
        "--unmapped-sequences",
        type=Path,
        required=True,
        help="Path to the unmapped_sequences.txt file.",
    )
    parser.add_argument(
        "--blast-results",
        type=Path,
        required=True,
        help="Path to blast results of unmapped_sequences.txt.",
    )
    parser.add_argument(
        "--config-file",
        type=Path,
        required=True,
        help="YAML file listing dataset configurations.",
    )
    parser.add_argument(
        "--blast-metadata",
        type=Path,
        required=True,
        help="Path to blast database metadata tsv file.",
    ),
    parser.add_argument(
        "--id-mapping",
        type=Path,
        required=False,
        help="Path to the ID mapping TSV file.",
    )
    parser.add_argument(
        "--excluded-sequences",
        type=Path,
        required=False,
        default=None,
        help=(
            "Path to the excluded_sequences.tsv file listing sequences kept out "
            "of the analysis for holding no usable base."
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Output file name.",
    )
    parser.add_argument(
        "--output-format",
        type=str,
        choices=["csv", "tsv", "json"],
        default="tsv",
        help="Output file format.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=False,
        default=None,
        help="ViralQC output directory (used to locate per-sample GFF and TBL files).",
    )
    args = parser.parse_args()

    # Load dataframes
    blast_metadata_df = load_blast_metadata(args.blast_metadata)
    formatted_dfs = format_dfs(args.files, args.config_file, blast_metadata_df)
    formatted_generic_dfs = format_dfs(
        args.generic_files, args.config_file, blast_metadata_df
    )
    unmapped_df = create_unmapped_df(
        args.unmapped_sequences, args.blast_results, blast_metadata_df
    )
    excluded_df = create_excluded_df(args.excluded_sequences)

    id_map = load_id_mapping(args.id_mapping) if args.id_mapping else {}
    reverse_id_map = {v: k for k, v in id_map.items()}

    # Resolve per-sample annotation directories
    gff_per_sample_dir = None
    tbl_per_sample_dir = None
    if args.output_dir:
        output_dir_abs = args.output_dir.resolve()
        gff_per_sample_dir = output_dir_abs / "gff_files" / "per_sample"
        tbl_per_sample_dir = output_dir_abs / "tbl_files" / "per_sample"

    # Organize sequences and results from Nextclade
    processed_seq_names = set()
    all_dfs = []

    for df in chain(formatted_dfs, formatted_generic_dfs):
        if "seqName" in df.columns:
            processed_seq_names.update(df["seqName"].astype(str).tolist())
        all_dfs.append(df)

    # Combine results into a single dataframe
    write_combined_df(
        iter(all_dfs),
        args.output,
        args.output_format,
        unmapped_df=unmapped_df,
        processed_seq_names=processed_seq_names,
        id_map=id_map,
        reverse_id_map=reverse_id_map,
        gff_per_sample_dir=gff_per_sample_dir,
        tbl_per_sample_dir=tbl_per_sample_dir,
        excluded_df=excluded_df,
    )
