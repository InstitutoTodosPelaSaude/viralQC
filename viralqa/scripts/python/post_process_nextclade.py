import argparse, re, csv, json
from pathlib import Path
from pandas import read_csv, concat, DataFrame
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
    "targetGene": str,
    "genomeQuality": str,
    "targetRegionsQuality": str,
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
    "coverage": "float64",
    "targetCdsCoverage": str,
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
        result[cds] = float(cov)
    return result


def get_target_regions_quality(
    cds_coverage: str, target_regions: list, target_regions_cov: float
) -> str:
    """
    Evaluate the quality of target regions. If any region has coverage
    lower than `target_regions_cov`, its status will be considered 'bad'.

    Args:
        cds_coverage: Value of the 'cdsCoverage' column from the Nextclade output.
        target_regions: List of target regions.
        target_regions_cov: Minimum required coverage for target regions.

    Returns:
        The status of the target regions.
    """
    cds_coverage = _parse_cds_cov(cds_coverage)
    if all(
        cds_coverage.get(region, 0) >= target_regions_cov for region in target_regions
    ):
        region_status = "good"
    else:
        region_status = "bad"

    return region_status


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
        virus_dataset = re.sub("\.nextclade.tsv", "", re.sub(".*\/", "", file))
        virus_info = config[virus_dataset]
        df = read_csv(file, sep="\t", header=0)
        df = format_sc2_clade(df, virus_dataset)
        df["virus"] = virus_info["virus_tag"]
        df["dataset"] = virus_info["dataset"]
        df["datasetVersion"] = virus_info["tag"]
        df["targetGene"] = ", ".join(virus_info["target_gene"])
        df["genomeQuality"] = df["qc.overallStatus"].apply(
            lambda x: StatusQuality[x].value
        )
        df["targetRegionsQuality"] = df["cdsCoverage"].apply(
            lambda cds_cov: get_target_regions_quality(
                cds_cov, virus_info["target_gene"], virus_info["target_gene_cov"]
            )
        )
        df["targetCdsCoverage"] = df["cdsCoverage"].apply(
            lambda cds_cov: get_target_regions_coverage(
                cds_cov, virus_info["target_gene"]
            )
        )
        dfs.append(df)

    return dfs


def create_unmapped_df(unmapped_sequences: Path) -> DataFrame:
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
    return df


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
    )

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
        "--files", nargs="+", help="List of Nextclade output .tsv files"
    )
    parser.add_argument(
        "--unmapped-sequences",
        type=Path,
        required=True,
        help="Path to the unmapped_sequences.txt file.",
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
    unmapped_df = create_unmapped_df(args.unmapped_sequences)
    formatted_dfs.append(unmapped_df)
    write_combined_df(formatted_dfs, args.output, args.output_format)
