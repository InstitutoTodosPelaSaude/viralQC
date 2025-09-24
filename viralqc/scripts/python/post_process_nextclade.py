import argparse, re, csv, os
from pathlib import Path
from pandas import read_csv, concat, DataFrame, notna
from pandas.errors import EmptyDataError
from yaml import safe_load
from enum import Enum


class StatusQuality(Enum):
    bad = "bad"
    mediocre = "median"
    good = "good"


TARGET_COLUMNS = {
    "seqName": str,
    "virus": str,
    "clade": str,
    "targetRegions": str,
    "targetGene": str,
    "genomeQuality": str,
    "targetRegionsQuality": str,
    "targetGeneQuality": str,
    "coverage": "float64",
    "cdsCoverage": str,
    "targetRegionsCoverage": str,
    "targetGeneCoverage": str,
    "qc.overallScore": "float64",
    "qc.overallStatus": str,
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
    "isReverseComplement": bool,
    "qc.missingData.missingDataThreshold": "float64",
    "qc.missingData.score": "float64",
    "qc.missingData.status": str,
    "qc.missingData.totalMissing": "Int64",
    "qc.snpClusters.score": "float64",
    "qc.snpClusters.status": str,
    "qc.snpClusters.totalSNPs": "Int64",
    "qc.frameShifts.totalFrameShifts": "Int64",
    "qc.frameShifts.totalFrameShiftsIgnored": "Int64",
    "qc.frameShifts.score": "float64",
    "qc.frameShifts.status": str,
    "qc.stopCodons.totalStopCodons": "Int64",
    "qc.stopCodons.score": "float64",
    "qc.stopCodons.status": str,
    "failedCdses": str,
    "warnings": str,
    "errors": str,
    "dataset": str,
    "datasetVersion": str,
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


def _parse_cds_cov(cds_list: str) -> list[dict[str, float]]:
    parts = cds_list.split(",")
    result = {}
    for p in parts:
        cds, cov = p.split(":")
        result[cds] = round(float(cov), 4)
    return result


def get_target_regions_quality(
    cds_coverage: str,
    genome_quality: str,
    target_regions: list,
    target_regions_cov: float,
) -> str:
    """
    Evaluate the quality of target regions. If any region has coverage
    lower than `target_regions_cov`, its status will be considered 'bad'.

    Args:
        cds_coverage: Value of the 'cdsCoverage' column from the Nextclade output.
        genome_quality: Quality of genome.
        target_regions: List of target regions.
        target_regions_cov: Minimum required coverage for target regions.

    Returns:
        The status of the target regions.
    """
    if genome_quality in ["good", ""]:
        return ""

    cds_coverage = _parse_cds_cov(cds_coverage)
    coverages = []
    for region in target_regions:
        coverages.append(float(cds_coverage.get(region, 0)))
    mean_coverage = sum(coverages) / len(coverages)

    if mean_coverage >= target_regions_cov:
        return "good"

    return "bad"


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
            virus_info = config[virus_dataset]
            df = format_sc2_clade(df, virus_dataset)
            df["virus"] = virus_info["virus_tag"]
            df["dataset"] = virus_info["dataset"]
            df["datasetVersion"] = virus_info["tag"]
            df["targetGene"] = virus_info["target_gene"]
            df["targetRegions"] = "-".join(virus_info["target_regions"])
            df["genomeQuality"] = df["qc.overallStatus"].apply(
                lambda x: StatusQuality[x].value if notna(x) else "missing"
            )
            df["targetRegionsQuality"] = df.apply(
                lambda row: (
                    get_target_regions_quality(
                        row["cdsCoverage"],
                        row["genomeQuality"],
                        virus_info["target_regions"],
                        virus_info["target_regions_cov"],
                    )
                    if notna(row["cdsCoverage"])
                    else ""
                ),
                axis=1,
            )
            df["targetRegionsCoverage"] = df["cdsCoverage"].apply(
                lambda cds_cov: (
                    get_target_regions_coverage(cds_cov, virus_info["target_regions"])
                    if notna(cds_cov)
                    else ""
                )
            )
            df["targetGeneQuality"] = df.apply(
                lambda row: (
                    get_target_regions_quality(
                        row["cdsCoverage"],
                        row["targetRegionsQuality"],
                        [virus_info["target_gene"]],
                        virus_info["target_gene_cov"],
                    )
                    if notna(row["cdsCoverage"])
                    else ""
                ),
                axis=1,
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
        dfs.append(df)

    return dfs


def _format_blast_virus_name(virus_name: str) -> str:
    formatted_virus_name = re.sub(".*_", "", virus_name)
    formatted_virus_name = re.sub("-", " ", formatted_virus_name)

    return formatted_virus_name


def create_unmapped_df(unmapped_sequences: Path, blast_results: Path) -> DataFrame:
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
        blast_df = read_csv(blast_results, sep="\t", header=None, names=blast_columns)
        blast_df_sub = blast_df[["seqName", "virus"]]

        merged = df.merge(
            blast_df_sub, on="seqName", how="left", suffixes=("_df1", "_df2")
        )
        merged["virus"] = merged["virus_df2"].combine_first(merged["virus_df1"])
        final_df = merged.drop(["virus_df1", "virus_df2"], axis=1)

        final_df = final_df.copy()
        final_df["virus"] = final_df["virus"].apply(_format_blast_virus_name)

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
    unmapped_df = create_unmapped_df(args.unmapped_sequences, args.blast_results)
    formatted_dfs.append(unmapped_df)
    write_combined_df(formatted_dfs, args.output, args.output_format)
