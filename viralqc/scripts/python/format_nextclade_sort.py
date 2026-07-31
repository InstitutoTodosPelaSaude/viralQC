from argparse import ArgumentParser
from pathlib import Path
from yaml import safe_load
from pandas import DataFrame, isna, read_csv, concat
from pandas.errors import EmptyDataError

# Schema of a nextclade sort results table, with the dtypes needed downstream.
NEXTCLADE_SORT_DTYPES = {
    "seqName": str,
    "dataset": str,
    "score": "Float64",
    "numHits": "Int64",
}
NEXTCLADE_SORT_COLUMNS = ["index", "seqName", "dataset", "score", "numHits"]


def read_nextclade_sort_tsv(tsv_file: Path) -> DataFrame:
    """
    Read a nextclade sort results table, tolerating an empty file.

    Handed a FASTA with no records, nextclade sort writes a zero-byte file with
    no header at all and still exits 0. That happens whenever every input
    sequence was excluded for holding no usable base, and read_csv raises
    EmptyDataError on such a file.

    Args:
        tsv_file: Path to the nextclade sort results TSV file.

    Returns:
        The parsed table, or an empty table carrying the expected columns.
    """
    try:
        return read_csv(tsv_file, sep="\t", dtype=NEXTCLADE_SORT_DTYPES)
    except EmptyDataError:
        return DataFrame(columns=NEXTCLADE_SORT_COLUMNS).astype(NEXTCLADE_SORT_DTYPES)


def map_datasets_to_local_paths(
    datasets_path: Path, config_file: Path
) -> dict[str, Path]:
    """
    Load dataset configuration from a YAML file and map each remote nextcladedataset name to a local path.

    Args:
        datasets_path: Base directory containing local dataset subdirectories.
        config_file: Path to the YAML configuration file listing nextclade datasets.

    Returns:
        A dictionary mapping remote dataset identifiers to their corresponding local directory paths.
    """
    with config_file.open("r") as f:
        config = safe_load(f)

    mapping: dict[str, Path] = {}
    for name, info in config["nextclade_data"].items():
        remote_name = info.get("dataset")
        if remote_name:
            mapping[remote_name] = datasets_path / name
    for name, info in config.get("github", {}).items():
        mapping[name] = datasets_path / name
    return mapping


def create_fasta_path(dataset: str, tsv_file: Path, seq_name: str) -> Path | None:
    """
    Construct the path to a FASTA file based on the TSV file location, dataset, and sequence filename.

    Args:
        dataset: Dataset identifier used in Nextclade output.
        tsv_file: Path to the Nextclade TSV results file.
        seq_name: Name of the sequences file as used in Nextclade.

    Returns:
        Path object pointing to the FASTA file location, or None if dataset is NaN.
    """
    if isna(dataset):
        return None
    return tsv_file.parent / dataset / seq_name


def format_nextclade_output(
    tsv_file: Path, tsv_external_file: Path, local_map: dict[str, Path]
) -> DataFrame:
    """
    Read and format Nextclade TSV output, attaching local dataset paths and FASTA file references.

    Args:
        tsv_file: Path to the Nextclade sort results TSV file.
        tsv_external_file: Path to the Nextclade sort results using external datasets TSV file.
        local_map: Mapping from remote dataset names to local paths.

    Returns:
        A DataFrame with added 'localDataset' and 'splittedFasta' columns.
    """
    df = read_nextclade_sort_tsv(tsv_file)
    df_external = read_nextclade_sort_tsv(tsv_external_file)
    df = concat([df, df_external], ignore_index=True)
    df = (
        df.sort_values(by="dataset", na_position="last")
        .drop_duplicates(subset="seqName", keep="first")
        .reset_index(drop=True)
    )
    df["localDataset"] = df["dataset"].map(local_map)
    df["splittedFasta"] = df["dataset"].apply(
        lambda dataset: create_fasta_path(dataset, tsv_file, "sequences.fa")
    )
    return df


def write_unmapped_sequences(formatted_df: DataFrame, output_dir: Path) -> None:
    """
    Identify sequences without a mapped dataset and write their names to a file.

    Args:
        formatted_df: Dataframe formatted by format_nextclade_output function
        output_dir: Output directory name.

    Returns:
        Nothing
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    unmapped_file = output_dir / "unmapped_sequences.txt"

    unmapped_seqs = formatted_df.loc[
        formatted_df["localDataset"].isna()
        | (formatted_df["localDataset"].astype(str).str.strip() == "")
        | (formatted_df["localDataset"].astype(str).str.lower() == "nan"),
        "seqName",
    ].drop_duplicates()
    if not unmapped_seqs.empty:
        unmapped_seqs.to_csv(unmapped_file, index=False, header=False)
    else:
        unmapped_file.write_text("")


if __name__ == "__main__":
    parser = ArgumentParser(
        description=(
            "Format Nextclade sort output and map entries to local dataset paths, "
            "then identify orphan sequences."
        )
    )
    parser.add_argument(
        "--nextclade-output",
        type=Path,
        required=True,
        help="Path to the Nextclade sort .tsv results file.",
    )
    parser.add_argument(
        "--nextclade-external-output",
        type=Path,
        required=True,
        help="Path to the Nextclade sort .tsv results file for external datasets.",
    ),
    parser.add_argument(
        "--config-file",
        type=Path,
        required=True,
        help="YAML file listing dataset configurations.",
    )
    parser.add_argument(
        "--local-datasets-path",
        type=Path,
        required=True,
        help="Base directory containing local datasets.",
    )
    parser.add_argument(
        "--output-path",
        type=Path,
        required=True,
        help="Directory to write output files.",
    )
    args = parser.parse_args()

    local_map = map_datasets_to_local_paths(args.local_datasets_path, args.config_file)
    formatted = format_nextclade_output(
        args.nextclade_output, args.nextclade_external_output, local_map
    )
    write_unmapped_sequences(formatted, args.output_path)

    formatted.dropna(subset=["dataset"]).to_csv(
        args.output_path / "datasets_selected.tsv", sep="\t", index=False
    )
