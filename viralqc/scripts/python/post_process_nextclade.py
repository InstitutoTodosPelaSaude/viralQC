import argparse, re, csv, os
from pathlib import Path
from pandas import read_csv, concat, DataFrame, notna, Series, NA
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
        df["clade"] = df["Nextclade_pango"]

    return df


def get_missing_data_quality(coverage: float) -> str:
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
    parts = cds_list.split(",")
    result = {}
    for p in parts:
        cds, cov = p.split(":")
        result[cds] = round(float(cov), 4)
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
    parts = cds_coverage.split(",")
    result = {}
    for p in parts:
        cds, cov = p.split(":")
        if float(cov) >= target_threshold_a:
            result[cds] = "A"
        elif float(cov) >= target_threshold_b:
            result[cds] = "B"
        elif float(cov) >= target_threshold_c:
            result[cds] = "C"
        elif float(cov) > 0:
            result[cds] = "D"

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

    cds_coverage = _parse_cds_cov(cds_coverage)
    cds_coverage = {k.strip(): v for k, v in cds_coverage.items()}
    coverages = []
    for region in target_regions:
        coverages.append(float(cds_coverage.get(region, 0)))
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
    df["targetRegionsCoverage"] = df["cdsCoverage"].apply(
        lambda cds_cov: (
            get_target_regions_coverage(cds_cov, virus_info["target_regions"])
            if notna(cds_cov)
            else ""
        )
    )
    df["targetGeneCoverage"] = df["cdsCoverage"].apply(
        lambda cds_cov: (
            get_target_regions_coverage(cds_cov, [virus_info["target_gene"]])
            if notna(cds_cov)
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
        missing_data_quality = get_missing_data_quality(row["coverage"])
        private_mutations_quality = get_private_mutations_quality(
            total=row["qc.privateMutations.total"],
            threshold=virus_info.get(
                "private_mutation_total_threshold",
                DEFAULT_PRIVATE_MUTATION_TOTAL_THRESHOLD,
            ),
        )
        mixed_sites_quality = get_qc_quality(row["qc.mixedSites.totalMixedSites"])
        snp_clusters_quality = get_qc_quality(row["qc.snpClusters.totalSNPs"])
        frameshifts_quality = get_qc_quality(row["qc.frameShifts.totalFrameShifts"])
        stop_codons_quality = get_qc_quality(row["qc.stopCodons.totalStopCodons"])

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
        if notna(row["cdsCoverage"]):
            target_regions_quality = get_target_regions_quality(
                cds_coverage=row["cdsCoverage"],
                genome_quality=genome_quality,
                target_regions=virus_info["target_regions"],
                target_threshold_a=COVERAGES_THRESHOLD["A"],
                target_threshold_b=COVERAGES_THRESHOLD["B"],
                target_threshold_c=COVERAGES_THRESHOLD["C"],
            )

            target_gene_quality = get_target_regions_quality(
                cds_coverage=row["cdsCoverage"],
                genome_quality=target_regions_quality,
                target_regions=[virus_info["target_gene"]],
                target_threshold_a=COVERAGES_THRESHOLD["A"],
                target_threshold_b=COVERAGES_THRESHOLD["B"],
                target_threshold_c=COVERAGES_THRESHOLD["C"],
            )

            cds_cov_quality = get_cds_cov_quality(
                cds_coverage=row["cdsCoverage"],
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
    return concat([df, qualities_df], axis=1)


def format_dfs(files: list[str], config_file: Path) -> list[DataFrame]:
    """
    Load and format nextclade outputs based on informations defined
    for each virus.

    Args:
        files: List of paths of nextclade outputs.
        config_file: Path to the YAML configuration file listing nextclade datasets.

    Returns:
        A list of formatted dataframes.
    """
    with config_file.open("r") as f:
        config = safe_load(f)
    dfs = []

    for file in files:
        try:
            df = read_csv(file, sep="\t", header=0)
        except EmptyDataError:
            df = DataFrame(columns=[TARGET_COLUMNS.keys()])

        if not df.empty:
            virus_dataset = re.sub("\.nextclade.tsv", "", re.sub(".*\/", "", file))
            virus_info = config["nextclade_data"].get(
                virus_dataset, config["github"].get(virus_dataset)
            )
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
        dfs.append(df)

    return dfs


def create_unmapped_df(
    unmapped_sequences: Path, blast_results: Path, blast_metadata: Path
) -> DataFrame:
    """
    Create a dataframe of unmapped sequences

    Args:
        unmapped_sequences: Path to unmapped_sequences.txt file
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
        return df
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
        names_metadata = [
            "virus",
            "segment",
            "virus_name",
            "virus_tax_id",
            "species_name",
            "species_tax_id",
            "dataset_with_version",
        ]

        blast_df = read_csv(blast_results, sep="\t", header=None, names=blast_columns)
        blast_metadata_df = read_csv(
            blast_metadata, sep="\t", header=None, names=names_metadata
        )

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

    return final_df


def write_combined_df(
    dfs: list[DataFrame], output_file: Path, output_format: str
) -> None:
    """
    Write a list of dataframes into a single file output.

    Args:
        dfs: A list of formatted dataframes.
        config_file: Path to output file
        output_format: format to write output (csv, tsv or json)

    Returns:
        Nothing
    """
    combined_df = concat(dfs, ignore_index=True)
    final_df = (
        combined_df[TARGET_COLUMNS.keys()]
        .astype(TARGET_COLUMNS)
        .sort_values(by=["virus"])
    ).round(4)
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

    formatted_dfs = format_dfs(args.files, args.config_file)
    unmapped_df = create_unmapped_df(
        args.unmapped_sequences, args.blast_results, args.blast_metadata
    )
    formatted_dfs.append(unmapped_df)
    write_combined_df(formatted_dfs, args.output, args.output_format)
