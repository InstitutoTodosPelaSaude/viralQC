import argparse, json
from pathlib import Path
from pandas import read_csv, concat, DataFrame
from glob import glob
from enum import Enum


class Separator(Enum):
    tsv = "\t"
    csv = ";"
    json = None


def read_gffs(files: list[str]) -> DataFrame:
    """
    Create a dataframe that represents the gff file of different viruses.

    Args:
        files: List of gff file paths.

    Returns:
        A dataframe that represents the GFF file.
    """
    column_names = [
        "seqname",
        "source",
        "feature",
        "start",
        "end",
        "score",
        "strand",
        "frame",
        "attribute",
    ]

    if not files:
        return DataFrame(columns=column_names)

    df = concat(
        (
            read_csv(
                f,
                delimiter="\t",
                comment="#",
                names=column_names,
                header=None,
            )
            for f in files
        ),
        ignore_index=True,
    )
    return df


def read_pp_nextclade(pp_results: Path, output_format: str) -> DataFrame:
    """
    Read results from post process nextclade independent of file format.

    Args:
        pp_results: Path to post process nextclade output file

    Returns:
        A dataframe that represents the post process nextclade output
    """
    sep = Separator[output_format].value
    if not sep:
        with open(pp_results) as f:
            js = json.load(f)
        df = DataFrame(js["data"]).set_index("index")
    else:
        df = read_csv(pp_results, sep=sep, header=0)
    return df


def check_target_regions(pp_results: DataFrame) -> dict:
    """
    Creates a dictionary mapping sequence names to target region names based on quality criteria.
    If the genomeQuality is "good", the target region will be "genome"; otherwise, it will take
    the value from the targetRegions or targetGene column.

    Briefly, the logic prioritizes regions with good quality status in the following order:
    genome > target regions > target genes.

    Args:
        pp_results: A DataFrame representing the post-processed Nextclade output.

    Returns:
        A dictionary with sequence names as keys and target region strings as values.
    """
    sequence_and_region = {}
    for _, row in pp_results.iterrows():
        if row["genomeQuality"] == "good":
            sequence_and_region[row["seqName"]] = "genome"
        elif row["targetRegionsQuality"] == "good":
            sequence_and_region[row["seqName"]] = row["targetRegions"]
        elif row["targetGeneQuality"] == "good":
            sequence_and_region[row["seqName"]] = row["targetGene"]
    return sequence_and_region


def get_regions(target_regions: dict, gff_info: DataFrame) -> dict:
    """
    Creates a dictionary with sequence name and target genomic positions. This function
    map the information presents on target regions dictionary with the information on a
    genomic features file (GFF) to create intervals of each region.

    Briefly, the logic look for start and end of gene features if the target region is
    a gene or a set of genes, or the start and end of the genome itsfel (feature = region)
    if the target region is the complete genome. 1 Is subtracted from start to produce a
    0-based bed.

    Args:
        target_regions: A dictionary with sequence name as key, and target regions as values.
        gff_info: A dataframe that represents the GFF file.
    Returns:
        A dictionary with sequence names as keys and a set containing the start and end positions
        of target regions.
    """
    sequences_intervals = {}
    seq_to_gff = {seq: df for seq, df in gff_info.groupby("seqname")}
    for seq, region in target_regions.items():
        df_seq = seq_to_gff.get(seq)
        if df_seq is None:
            continue

        if region == "genome":
            region_rows = df_seq[df_seq["feature"] == "region"]
            if not region_rows.empty:
                sequences_intervals[seq] = (
                    region_rows["start"].values[0] - 1,
                    region_rows["end"].values[0],
                )
        else:
            genes = region.split("-")
            gene_rows = df_seq[df_seq["feature"] == "gene"]
            mask = gene_rows["attribute"].str.contains(
                "|".join(f"gene_name={g}" for g in genes)
            )
            gene_rows = gene_rows[mask]

            if not gene_rows.empty:
                sequences_intervals[seq] = (
                    gene_rows["start"].min() - 1,
                    gene_rows["end"].max(),
                )
    return sequences_intervals


def write_bed(sequences_intervals: dict, output_file: Path) -> None:
    """
    Writes the target regions intervals to a bed file.

    Args:
        sequences_intervals: A dictionary with sequence names as keys and a set
        containing the start and end positions of target regions.
        output_file: Output file path.
    """
    with output_file.open("w") as f:
        for seqname, (start, end) in sequences_intervals.items():
            f.write(f"{seqname}\t{int(start)}\t{int(end)}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generates a bed file based on regions considered with good quality by nextclade."
    )
    parser.add_argument(
        "--pp-results",
        type=Path,
        required=True,
        help="Path to the post process nextclade task output file.",
    )
    parser.add_argument(
        "--output-format",
        type=str,
        choices=["csv", "tsv", "json"],
        default="tsv",
        help="Post process nextclade output file format.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Output file name.",
    )
    args = parser.parse_args()

    files = glob(f"{args.pp_results.parent}/*.gff")
    gff = read_gffs(files)
    pp_nextclade = read_pp_nextclade(args.pp_results, args.output_format)
    target_regions = check_target_regions(pp_nextclade)
    target_regions_intervals = get_regions(target_regions, gff)
    write_bed(target_regions_intervals, args.output)
