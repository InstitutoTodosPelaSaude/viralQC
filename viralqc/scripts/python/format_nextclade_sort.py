from argparse import ArgumentParser
from pathlib import Path
from yaml import safe_load
from pandas import DataFrame, isna, read_csv


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
    for name, info in config.items():
        remote_name = info.get("dataset")
        if remote_name:
            mapping[remote_name] = datasets_path / name
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


def format_nextclade_output(tsv_file: Path, local_map: dict[str, Path]) -> DataFrame:
    """
    Read and format Nextclade TSV output, attaching local dataset paths and FASTA file references.

    Args:
        tsv_file: Path to the Nextclade sort results TSV file.
        local_map: Mapping from remote dataset names to local paths.

    Returns:
        A DataFrame with added 'localDataset' and 'splittedFasta' columns.
    """
    df = read_csv(
        tsv_file,
        sep="\t",
        index_col=0,
        dtype={"seqName": str, "dataset": str, "score": "Float64", "numHits": "Int64"},
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

    unmapped_seqs = formatted_df.loc[formatted_df["localDataset"].isna(), "seqName"]
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
    formatted = format_nextclade_output(args.nextclade_output, local_map)
    write_unmapped_sequences(formatted, args.output_path)

    formatted.dropna(subset=["dataset"]).to_csv(
        args.output_path / "datasets_selected.tsv", sep="\t", index=False
    )
