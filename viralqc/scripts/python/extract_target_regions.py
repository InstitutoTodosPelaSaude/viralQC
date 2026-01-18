import argparse, json, gc
from pathlib import Path
from typing import Generator, Iterator
from pandas import read_csv, concat, DataFrame
from glob import glob
from enum import Enum


class Separator(Enum):
    tsv = "\t"
    csv = ";"
    json = None


def load_id_mapping(mapping_path: Path) -> dict:
    """
    Load the ID mapping TSV file. Returns a dict mapping Original Header -> Sanitized ID.
    """
    mapping = {}
    try:
        with open(mapping_path, "r") as f:
            reader = read_csv(f, sep="\t", dtype=str)
            return dict(zip(reader["original_header"], reader["id"]))
    except Exception:
        return {}


# Columns required from post-process nextclade file with their dtypes
REQUIRED_COLUMNS = {
    "seqName": str,
    "genomeQuality": str,
    "targetRegionsQuality": str,
    "targetGeneQuality": str,
    "targetRegions": str,
    "targetGene": str,
}

DEFAULT_CHUNK_SIZE = 10000


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
                dtype={"seqname": str},
            )
            for f in files
        ),
        ignore_index=True,
    )
    return df


def read_pp_nextclade_chunks(
    pp_results: Path, output_format: str, chunk_size: int = DEFAULT_CHUNK_SIZE
) -> Generator[DataFrame, None, None]:
    """
    Read results from post process nextclade in chunks for memory efficiency.

    Args:
        pp_results: Path to post process nextclade output file
        output_format: Format of the input file (csv, tsv, or json)
        chunk_size: Number of rows to read per chunk

    Yields:
        DataFrame chunks from the post process nextclade output
    """
    sep = Separator[output_format].value

    if not sep:
        with open(pp_results) as f:
            js = json.load(f)

        df = DataFrame(js["data"])
        del js
        gc.collect()

        available_cols = [col for col in REQUIRED_COLUMNS.keys() if col in df.columns]
        df = df[available_cols]

        for start in range(0, len(df), chunk_size):
            chunk = df.iloc[start : start + chunk_size].copy()
            yield chunk
            del chunk
            gc.collect()

        del df
        gc.collect()
    else:
        chunks = read_csv(
            pp_results,
            sep=sep,
            header=0,
            usecols=lambda col: col in REQUIRED_COLUMNS,
            dtype=REQUIRED_COLUMNS,
            chunksize=chunk_size,
        )

        for chunk in chunks:
            yield chunk
            del chunk
            gc.collect()


def get_region_interval(
    seq: str, region: str, seq_to_gff: dict[str, DataFrame]
) -> tuple | None:
    """
    Get the genomic interval for a single sequence and region.

    Args:
        seq: Sequence name
        region: Target region string
        seq_to_gff: Dictionary mapping sequence names to GFF dataframes

    Returns:
        Tuple of (start, end, region_name) or None if not found
    """
    if " " in seq:
        seq_norm = seq.split()[0]
    else:
        seq_norm = seq

    df_seq = seq_to_gff.get(seq_norm)
    if df_seq is None:
        return None

    if isinstance(region, float):
        return None

    if region == "genome":
        region_rows = df_seq[df_seq["feature"] == "region"]
        if not region_rows.empty:
            return (
                region_rows["start"].values[0] - 1,
                region_rows["end"].values[0],
                "genome",
            )
    else:
        genes = region.split("|")
        gene_rows = df_seq[df_seq["feature"] == "gene"]
        mask = gene_rows["attribute"].str.contains(
            "|".join(f"gene_name={g}" for g in genes)
        )
        gene_rows = gene_rows[mask]

        if not gene_rows.empty:
            region_name = ",".join(genes)
            return (
                gene_rows["start"].min() - 1,
                gene_rows["end"].max(),
                region_name,
            )

    return None


def get_sanitized_id(seq_name: str, id_map: dict) -> str:
    """
    Get the sanitized ID for a sequence name.
    If not in map, return original (fallback).
    """
    return id_map.get(seq_name, seq_name)


def process_and_write_bed(
    chunks_iter: Iterator[DataFrame],
    gff_info: DataFrame,
    output_file: Path,
    id_map: dict,
) -> None:
    """
    Process chunks and write BED file incrementally without accumulating data in memory.

    Args:
        chunks_iter: Iterator yielding DataFrame chunks
        gff_info: A dataframe that represents the GFF file
        output_file: Output file path
        id_map: Dictionary mapping Original Header -> Sanitized ID
    """
    seq_to_gff = {seq: df for seq, df in gff_info.groupby("seqname")}
    written_sequences = set()

    with output_file.open("w") as f:
        for chunk in chunks_iter:
            for _, row in chunk.iterrows():
                seq_name = row["seqName"]

                if seq_name in written_sequences:
                    continue

                sanitized_id = (
                    get_sanitized_id(seq_name, id_map) if id_map else seq_name
                )

                genome_quality = row.get("genomeQuality", "")
                target_regions_quality = row.get("targetRegionsQuality", "")
                target_gene_quality = row.get("targetGeneQuality", "")

                region = None
                if genome_quality in ["A", "B"]:
                    region = "genome"
                elif target_regions_quality in ["A", "B"]:
                    region = row.get("targetRegions", "")
                elif target_gene_quality in ["A", "B"]:
                    region = row.get("targetGene", "")

                if region is None:
                    continue

                interval = get_region_interval(sanitized_id, region, seq_to_gff)

                if interval is not None:
                    start, end, region_name = interval
                    f.write(f"{seq_name}\t{int(start)}\t{int(end)}\t{region_name}\n")
                    written_sequences.add(seq_name)

            del chunk
            gc.collect()

    del seq_to_gff, written_sequences
    gc.collect()


# Keep original functions for backwards compatibility
def read_pp_nextclade(pp_results: Path, output_format: str) -> DataFrame:
    """
    Read results from post process nextclade independent of file format.

    Args:
        pp_results: Path to post process nextclade output file
        output_format: Format of the input file

    Returns:
        A dataframe that represents the post process nextclade output
    """
    sep = Separator[output_format].value
    if not sep:
        with open(pp_results) as f:
            js = json.load(f)
        df = DataFrame(js["data"]).set_index("index")
        available_cols = [col for col in REQUIRED_COLUMNS.keys() if col in df.columns]
        df = df[available_cols]
    else:
        df = read_csv(
            pp_results,
            sep=sep,
            header=0,
            usecols=lambda col: col in REQUIRED_COLUMNS,
            dtype=REQUIRED_COLUMNS,
        )
    return df


def check_target_regions(pp_results: DataFrame) -> dict:
    """
    Creates a dictionary mapping sequence names to target region names based on quality criteria.
    """
    sequence_and_region = {}
    for _, row in pp_results.iterrows():
        if row["genomeQuality"] in ["A", "B"]:
            sequence_and_region[row["seqName"]] = "genome"
        elif row["targetRegionsQuality"] in ["A", "B"]:
            sequence_and_region[row["seqName"]] = row["targetRegions"]
        elif row["targetGeneQuality"] in ["A", "B"]:
            sequence_and_region[row["seqName"]] = row["targetGene"]
    return sequence_and_region


def get_regions(target_regions: dict, gff_info: DataFrame) -> dict:
    """
    Creates a dictionary with sequence name and target genomic positions.
    """
    sequences_intervals = {}
    seq_to_gff = {seq: df for seq, df in gff_info.groupby("seqname")}
    for seq, region in target_regions.items():
        interval = get_region_interval(seq, region, seq_to_gff)
        if interval is not None:
            sequences_intervals[seq] = interval
    return sequences_intervals


def write_bed(sequences_intervals: dict, output_file: Path) -> None:
    """
    Writes the target regions intervals to a bed file with region names.
    """
    with output_file.open("w") as f:
        for seqname, (start, end, region_name) in sequences_intervals.items():
            f.write(f"{seqname}\t{int(start)}\t{int(end)}\t{region_name}\n")


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
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=DEFAULT_CHUNK_SIZE,
        help="Number of rows to process at a time (default: 10000).",
    )
    parser.add_argument(
        "--id-mapping",
        type=Path,
        required=False,
        help="Path to the ID mapping TSV file.",
    )
    args = parser.parse_args()

    files = glob(f"{args.pp_results.parent}/gff_files/*.gff")
    gff = read_gffs(files)

    id_map = load_id_mapping(args.id_mapping) if args.id_mapping else {}

    pp_chunks = read_pp_nextclade_chunks(
        args.pp_results, args.output_format, args.chunk_size
    )
    process_and_write_bed(pp_chunks, gff, args.output, id_map)

    del gff
    gc.collect()
