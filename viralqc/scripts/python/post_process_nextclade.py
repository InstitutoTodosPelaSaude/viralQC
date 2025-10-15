import argparse, re, csv, os
from pathlib import Path
from pandas import read_csv, concat, DataFrame, notna
from pandas.errors import EmptyDataError
from yaml import safe_load


target_threshold_cOLUMNS = {
    "seqName": str,
    "virus": str,
    "segment": str,
    "ncbi_id": str,
    "clade": str,
    "targetRegions": str,
    "targetGene": str,
    "genomeQuality": str,
    "targetRegionsQuality": str,
    "targetGeneQuality": str,
    "cdsCoverageQuality": str,
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
    "qc.missingData.score": "float64",
    "qc.missingData.status": str,
    "qc.snpClusters.score": "float64",
    "qc.snpClusters.status": str,
    "qc.frameShifts.score": "float64",
    "qc.frameShifts.status": str,
    "qc.stopCodons.score": "float64",
    "qc.stopCodons.status": str,
    "dataset": str,
    "datasetVersion": str,
}

# DEFAULT COVERAGE AND QUALITY SCORES
GENOME_COV = 0.7
GEONOME_SCORES = {
    "A": 50,
    "B": 75,
    "C": 100,
}

GENE_COVERAGES = {
    "A": 0.95,
    "B": 0.75,
    "C": 0.5,
}

REGION_COVERAGES = {
    "A": 0.95,
    "B": 0.75,
    "C": 0.50,
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


def get_genome_quality(
    qc_overall_score: float,
    genome_coverage: float,
    genome_coverage_threshold: float,
    genome_score_threshold_a: float,
    genome_score_threshold_b: float,
    genome_score_threshold_c: float,
) -> str:
    """
    Evaluate the quality of genome and classify it as categories based
    on coverage and score,

    Args:
        qc_overall_score: Value of the 'qc.OverallScore' column from the Nextclade output.
        genome_coverage: Value of the 'coverage' column from the Nextclade output.
        genome_coverage_threshold: Minimum required coverage for consider a genome as "A".
        genome_score_threshold_a: Minimum required coverage for consider a target regions as "A".
        genome_score_threshold_b: Minimum required coverage for consider a target regions as "B".
        genome_score_threshold_c: Minimum required coverage for consider a target regions as "C".

    Returns:
        The quality of genome
    """
    if qc_overall_score is None:
        return None
    elif (
        genome_coverage >= genome_coverage_threshold
        and qc_overall_score <= genome_score_threshold_a
    ):
        return "A"
    elif qc_overall_score <= genome_score_threshold_b:
        return "B"
    elif qc_overall_score <= genome_score_threshold_c:
        return "C"
    else:
        return "D"


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
    if genome_quality in ["A", ""]:
        return ""

    cds_coverage = _parse_cds_cov(cds_coverage)
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
    target_threshold_cds_coverage = [
        f"{region}: {cds_coverage.get(region,0)}" for region in target_regions
    ]

    return ", ".join(target_threshold_cds_coverage)


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
            df = DataFrame(columns=[target_threshold_cOLUMNS.keys()])

        if not df.empty:
            virus_dataset = re.sub("\.nextclade.tsv", "", re.sub(".*\/", "", file))
            virus_info = config[virus_dataset]
            df = format_sc2_clade(df, virus_dataset)
            df["virus"] = virus_info["virus_name"]
            df["segment"] = virus_info["segment"]
            df["ncbi_id"] = virus_info["ncbi_id"]
            df["dataset"] = virus_info["dataset"]
            df["datasetVersion"] = virus_info["tag"]
            df["targetGene"] = virus_info["target_gene"]
            df["targetRegions"] = "|".join(virus_info["target_regions"])
            # df["genomeQuality"] = df["qc.overallStatus"].apply(
            #     lambda x: StatusQuality[x].value if notna(x) else "missing"
            # )
            df["genomeQuality"] = df.apply(
                lambda row: (
                    get_genome_quality(
                        qc_overall_score=row["qc.overallScore"],
                        genome_coverage=row["coverage"],
                        genome_coverage_threshold=virus_info.get(
                            "genome_cov", GENOME_COV
                        ),
                        genome_score_threshold_a=virus_info.get(
                            "genome_score_threshold", GEONOME_SCORES
                        ).get("A"),
                        genome_score_threshold_b=virus_info.get(
                            "genome_score_threshold", GEONOME_SCORES
                        ).get("B"),
                        genome_score_threshold_c=virus_info.get(
                            "genome_score_threshold", GEONOME_SCORES
                        ).get("C"),
                    )
                ),
                axis=1,
            )
            df["targetRegionsQuality"] = df.apply(
                lambda row: (
                    get_target_regions_quality(
                        cds_coverage=row["cdsCoverage"],
                        genome_quality=row["genomeQuality"],
                        target_regions=virus_info["target_regions"],
                        target_threshold_a=virus_info.get(
                            "target_regions_cov", REGION_COVERAGES
                        ).get("A"),
                        target_threshold_b=virus_info.get(
                            "target_regions_cov", REGION_COVERAGES
                        ).get("B"),
                        target_threshold_c=virus_info.get(
                            "target_regions_cov", REGION_COVERAGES
                        ).get("C"),
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
                        cds_coverage=row["cdsCoverage"],
                        genome_quality=row["targetRegionsQuality"],
                        target_regions=[virus_info["target_gene"]],
                        target_threshold_a=virus_info.get(
                            "target_regions_cov", GENE_COVERAGES
                        ).get("A"),
                        target_threshold_b=virus_info.get(
                            "target_regions_cov", GENE_COVERAGES
                        ).get("B"),
                        target_threshold_c=virus_info.get(
                            "target_regions_cov", GENE_COVERAGES
                        ).get("C"),
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
            df["cdsCoverageQuality"] = df.apply(
                lambda row: (
                    get_cds_cov_quality(
                        cds_coverage=row["cdsCoverage"],
                        target_threshold_a=virus_info.get(
                            "target_regions_cov", REGION_COVERAGES
                        ).get("A"),
                        target_threshold_b=virus_info.get(
                            "target_regions_cov", REGION_COVERAGES
                        ).get("B"),
                        target_threshold_c=virus_info.get(
                            "target_regions_cov", REGION_COVERAGES
                        ).get("C"),
                    )
                    if notna(row["cdsCoverage"])
                    else ""
                ),
                axis=1,
            )
        dfs.append(df)

    return dfs


def _format_blast_virus_name(virus_name: str) -> str:
    formatted_virus_name = re.sub(".*_", "", virus_name)
    formatted_virus_name = re.sub("-", " ", formatted_virus_name)

    return formatted_virus_name


def _get_blast_virus_id(virus_name: str) -> str:
    parts = virus_name.split("_")
    virus_id = parts[0] + "_" + parts[1]

    return virus_id


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

    for col in target_threshold_cOLUMNS.keys():
        if col not in df.columns:
            if target_threshold_cOLUMNS[col] == str:
                df[col] = ""
            elif target_threshold_cOLUMNS[col] == "float64":
                df[col] = None
            elif target_threshold_cOLUMNS[col] == "Int64":
                df[col] = None
            elif target_threshold_cOLUMNS[col] == bool:
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
        final_df["ncbi_id"] = final_df["virus"].apply(_get_blast_virus_id)
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
        combined_df[target_threshold_cOLUMNS.keys()]
        .astype(target_threshold_cOLUMNS)
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
