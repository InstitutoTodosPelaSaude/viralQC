import argparse, re, csv, os
from pathlib import Path
from pandas import read_csv, concat, DataFrame, notna, Series, NA, to_numeric
from numpy import nan
from pandas.errors import EmptyDataError
from yaml import safe_load


TARGET_COLUMNS = {
    "seqName": str,
    "virus": str,
    "virus_tax_id": "Int64",
    "virus_species": str,
    "virus_species_tax_id": "Int64",
    "segment": str,
    "ncbi_id": str,
    "clade": str,
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
    "dataset": str,
    "datasetVersion": str,
}


DEFAULT_PRIVATE_MUTATION_TOTAL_THRESHOLD = 10
COVERAGES_THRESHOLD = {
    "A": 0.95,
    "B": 0.75,
    "C": 0.5,
}


def load_blast_metadata(metadata_path: Path) -> DataFrame:
    """
    Load the BLAST metadata TSV file.

    Args:
        metadata_path: Path to the metadata file.

    Returns:
        Dataframe containing the BLAST metadata.
    """
    names_metadata = [
        "virus",
        "segment",
        "virus_name",
        "virus_tax_id",
        "species_name",
        "species_tax_id",
        "dataset_with_version",
    ]
    try:
        df = read_csv(metadata_path, sep="\t", header=None, names=names_metadata)
        return df.loc[:, ~df.columns.duplicated()]
    except Exception:
        return DataFrame(columns=names_metadata)


def format_sc2_clade(df: DataFrame, dataset_name: str) -> DataFrame:
    """
    For SARS-CoV-2 datasets, replaces 'clade' with 'Nextclade_pango'.

    Args:
        df: Dataframe of nextclade results.
        dataset_name: Name of dataset.

    Returns:
        For SARS-CoV-2 datasets returns a dataframe with values from
        Nextclade_pango column into clade column.
    """
    if dataset_name.startswith("sarscov2"):
        df = df.copy()
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
        The quality of genome
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


def _parse_cds_cov(cds_list: str) -> list[dict[str, float]]:
    """
    Parse the cdsCoverage string into a dictionary.

    Args:
        cds_list: String containing CDS coverage data.

    Returns:
        Dictionary mapping gene names to coverage values.
    """
    if not isinstance(cds_list, str):
        return {}
    parts = cds_list.split(",")
    result = {}
    for p in parts:
        if ":" in p:
            cds, cov = p.split(":")
            try:
                result[cds] = round(float(cov), 4)
            except ValueError:
                continue
    return result


def get_cds_cov_quality(
    cds_coverage: str,
    target_threshold_a: float,
    target_threshold_b: float,
    target_threshold_c: float,
) -> list[dict[str, str]]:
    """
    Categorize the cds regions based on coverage thresholds.

    Args:
        cds_coverage: Value of the 'cdsCoverage' column from the Nextclade output.
        target_threshold_a: Minimum required coverage for consider a target regions as "A".
        target_threshold_b: Minimum required coverage for consider a target regions as "B".
        target_threshold_c: Minimum required coverage for consider a target regions as "C".

    Returns:
        The status of the target regions.
    """
    if not isinstance(cds_coverage, str):
        return ""

    parts = cds_coverage.split(",")
    result = {}
    for p in parts:
        if ":" in p:
            cds, cov = p.split(":")
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
            except ValueError:
                continue

    return ", ".join(f"{cds}: {coverage}" for cds, coverage in result.items())


def get_target_regions_quality(
    cds_coverage: str,
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
        cds_coverage: Value of the 'cdsCoverage' column from the Nextclade output.
        genome_quality: Quality of genome.
        target_regions: List of target regions.
        target_threshold_a: Minimum required coverage for consider a target regions as "A".
        target_threshold_b: Minimum required coverage for consider a target regions as "B".
        target_threshold_c: Minimum required coverage for consider a target regions as "C".

    Returns:
        The status of the target regions.
    """
    if genome_quality in ["A", "B", ""]:
        return ""

    if not target_regions:
        return ""

    cds_coverage = _parse_cds_cov(cds_coverage)
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


def get_target_regions_coverage(cds_coverage: str, target_regions: list[str]) -> str:
    """
    Extract the coverage of specific genomic regions.

    Args:
        cds_coverage: Value of the 'cdsCoverage' column from the Nextclade output.
        target_regions: List of target regions.

    Returns:
        A string with region and coverage.
    """
    cds_coverage = _parse_cds_cov(cds_coverage)
    target_cds_coverage = [
        f"{region}: {cds_coverage.get(region,0)}" for region in target_regions
    ]

    return ", ".join(target_cds_coverage)


def add_coverages(df: DataFrame, virus_info: dict) -> DataFrame:
    """
    Add 'targetRegionsCoverage', 'targetGeneCoverage' and format
    'cdsCoverage' column to results datafarame.

    Args:
        df: Dataframe of nextclade results.
        virus_info: Dictionary with specific virus configuration

    Returns:
        The dataframe with the new columns.
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

    df["cdsCoverage"] = df["cdsCoverage"].apply(_parse_cds_cov)
    df["cdsCoverage"] = df["cdsCoverage"].apply(
        lambda d: ", ".join(f"{cds}: {coverage}" for cds, coverage in d.items())
    )
    return df


def add_qualities(df: DataFrame, virus_info: dict) -> DataFrame:
    """
    Compute all quality metrics into a single apply.

    Args:
        df: Dataframe of nextclade results.
        virus_info: Dictionary with specific virus configuration

    Returns:
        The dataframe with the new quality columns.
    """

    def compute_all_qualities(row):
        # --- Metrics qualities ---
        missing_data_quality = get_missing_data_quality(row.get("coverage", nan))

        private_mutations_total = row.get("qc.privateMutations.total", nan)
        private_mutations_quality = get_private_mutations_quality(
            total=private_mutations_total,
            threshold=virus_info.get(
                "private_mutation_total_threshold",
                DEFAULT_PRIVATE_MUTATION_TOTAL_THRESHOLD,
            ),
        )

        mixed_sites_quality = get_qc_quality(
            row.get("qc.mixedSites.totalMixedSites", nan)
        )
        snp_clusters_quality = get_qc_quality(row.get("qc.snpClusters.totalSNPs", nan))
        frameshifts_quality = get_qc_quality(
            row.get("qc.frameShifts.totalFrameShifts", nan)
        )
        stop_codons_quality = get_qc_quality(
            row.get("qc.stopCodons.totalStopCodons", nan)
        )

        # --- Genome quality ---
        genome_score, genome_quality = get_genome_quality(
            [
                missing_data_quality,
                mixed_sites_quality,
                private_mutations_quality,
                snp_clusters_quality,
                frameshifts_quality,
                stop_codons_quality,
            ]
        )

        # --- Target qualities ---
        cds_coverage = row.get("cdsCoverage", nan)
        if notna(cds_coverage) and cds_coverage != "":
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
                target_threshold_a=virus_info.get(
                    "target_regions_cov", COVERAGES_THRESHOLD
                )["A"],
                target_threshold_b=virus_info.get(
                    "target_regions_cov", COVERAGES_THRESHOLD
                )["B"],
                target_threshold_c=virus_info.get(
                    "target_regions_cov", COVERAGES_THRESHOLD
                )["C"],
            )
        else:
            target_regions_quality = ""
            target_gene_quality = ""
            cds_cov_quality = ""

        return Series(
            {
                "missingDataQuality": missing_data_quality,
                "privateMutationsQuality": private_mutations_quality,
                "mixedSitesQuality": mixed_sites_quality,
                "snpClustersQuality": snp_clusters_quality,
                "frameShiftsQuality": frameshifts_quality,
                "stopCodonsQuality": stop_codons_quality,
                "genomeQualityScore": genome_score,
                "genomeQuality": genome_quality,
                "targetRegionsQuality": target_regions_quality,
                "targetGeneQuality": target_gene_quality,
                "cdsCoverageQuality": cds_cov_quality,
            }
        )

    qualities_df = df.apply(compute_all_qualities, axis=1)

    # Drop columns that are about to be added to avoid duplicates
    df = df.drop(columns=qualities_df.columns, errors="ignore")

    return concat([df, qualities_df], axis=1)


def format_dfs(
    files: list[str], config_file: Path, blast_metadata_df: DataFrame = None
) -> list[DataFrame]:
    """
    Load and format nextclade outputs based on informations defined
    for each virus.

    Args:
        files: List of paths of nextclade outputs.
        config_file: Path to the YAML configuration file listing nextclade datasets.
        blast_metadata_df: Dataframe with BLAST metadata (optional).

    Returns:
        A list of formatted dataframes.
    """
    with config_file.open("r") as f:
        config = safe_load(f)
    dfs = []

    for file in files:
        try:
            df = read_csv(file, sep="\t", header=0)
            df = df.loc[:, ~df.columns.duplicated()]
        except EmptyDataError:
            df = DataFrame(columns=list(TARGET_COLUMNS.keys()))

        if not df.empty:
            virus_dataset = re.sub("\.nextclade.tsv", "", re.sub(".*\/", "", file))
            virus_dataset = re.sub("\.generic", "", virus_dataset)

            virus_info = config["nextclade_data"].get(
                virus_dataset, config["github"].get(virus_dataset)
            )

            if virus_info:
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
                df = add_coverages(df, virus_info)
                df = add_qualities(df, virus_info)
            else:
                # Nextclade generic run
                # virus_dataset is the accession (e.g., AC_000006.1)
                df["ncbi_id"] = virus_dataset

                # Enrich with BLAST metadata if available
                if blast_metadata_df is not None:
                    # Filter metadata for this accession
                    meta = blast_metadata_df[
                        blast_metadata_df["virus"] == virus_dataset
                    ]
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

                # Explicitly set columns to nan/empty for generic runs
                df["clade"] = nan
                df["qc.overallScore"] = nan
                df["qc.overallStatus"] = nan

                # Add empty columns for missing info
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

                # For generic runs, we DO NOT calculate qualities as we don't have thresholds
                # But we DO want to format cdsCoverage if it exists

                mock_virus_info = {
                    "target_regions": [],
                    "target_gene": "",
                }
                df = add_coverages(df, mock_virus_info)

                # Explicitly set quality columns to empty/None
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

        # Ensure no duplicates before appending
        df = df.loc[:, ~df.columns.duplicated()]
        dfs.append(df)

    return dfs


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
    with open(unmapped_sequences, "r") as f:
        data = [(line.strip(), "Unclassified") for line in f]
    df = DataFrame(data, columns=["seqName", "virus"])

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

        merged = df.merge(blast_df, on="seqName", how="left", suffixes=("_df1", "_df2"))
        merged["virus"] = merged["virus_df2"].fillna(merged["virus_df1"])

        final_df = merged.drop(columns=["virus_df1", "virus_df2"])
        final_df = final_df.assign(
            ncbi_id=final_df["virus"],
            virus=final_df["virus_name"].fillna("Unclassified").astype(str),
            virus_tax_id=final_df["virus_tax_id_df2"].astype("Int64"),
            virus_species=final_df["species_name"].fillna("Unclassified").astype(str),
            virus_species_tax_id=final_df["species_tax_id"].astype("Int64"),
            segment=final_df["segment_df2"].fillna("Unsegmented").astype(str),
        )
        final_df[["dataset", "datasetVersion"]] = final_df[
            "dataset_with_version"
        ].str.split("_", n=1, expand=True)

    return final_df.loc[:, ~final_df.columns.duplicated()]


def write_combined_df(
    dfs: list[DataFrame], output_file: Path, output_format: str
) -> None:
    """
    Write a list of dataframes into a single file output.

    Args:
        dfs: A list of formatted dataframes.
        output_file: Path to output file
        output_format: format to write output (csv, tsv or json)

    Returns:
        Nothing
    """
    # Ensure all dfs have unique columns before concat
    dfs = [df.loc[:, ~df.columns.duplicated()] for df in dfs]

    combined_df = concat(dfs, ignore_index=True)

    # Ensure combined_df has unique columns
    combined_df = combined_df.loc[:, ~combined_df.columns.duplicated()]

    # Sanitize columns to ensure they can be cast safely
    for col, dtype in TARGET_COLUMNS.items():
        if col in combined_df.columns:
            if dtype == "Int64" or dtype == "float64":
                combined_df[col] = to_numeric(combined_df[col], errors="coerce").astype(
                    dtype
                )
            elif dtype == str:
                # Use pandas nullable string type to preserve NA/None
                combined_df[col] = combined_df[col].astype("string")

    # Select columns and sort
    final_df = (
        combined_df[list(TARGET_COLUMNS.keys())].sort_values(by=["virus"])
    ).round(4)

    # Convert whitespace strings to NaN, which will be handled correctly by output formats
    # (empty string for CSV/TSV, null for JSON)
    final_df = final_df.replace(r"^\s*$", nan, regex=True)

    if output_format == "tsv":
        final_df.to_csv(output_file, sep="\t", index=False, header=True)
    if output_format == "csv":
        final_df.to_csv(
            output_file, sep=";", index=False, header=True, quoting=csv.QUOTE_NONNUMERIC
        )
    if output_format == "json":
        json_content = final_df.to_json(orient="table", indent=4)
        json_content = json_content.replace("\\/", "/")
        with open(output_file, "w", encoding="utf-8") as f:
            f.write(json_content)


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
        help="Output file name.",
    )
    args = parser.parse_args()

    # Load BLAST metadata once
    blast_metadata_df = load_blast_metadata(args.blast_metadata)

    formatted_dfs = format_dfs(args.files, args.config_file, blast_metadata_df)
    formatted_generic_dfs = format_dfs(
        args.generic_files, args.config_file, blast_metadata_df
    )

    # Combine all Nextclade results (standard + generic)
    all_nextclade_dfs = formatted_dfs + formatted_generic_dfs

    # Collect all sequence names that have been processed by Nextclade
    processed_seq_names = set()
    for df in all_nextclade_dfs:
        if "seqName" in df.columns:
            processed_seq_names.update(df["seqName"].astype(str).tolist())

    unmapped_df = create_unmapped_df(
        args.unmapped_sequences, args.blast_results, blast_metadata_df
    )

    # Filter unmapped_df to exclude sequences that were actually processed
    if not unmapped_df.empty and "seqName" in unmapped_df.columns:
        unmapped_df = unmapped_df[
            ~unmapped_df["seqName"].astype(str).isin(processed_seq_names)
        ]

    all_nextclade_dfs.append(unmapped_df)
    write_combined_df(all_nextclade_dfs, args.output, args.output_format)
