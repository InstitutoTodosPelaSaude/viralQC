"""Create standard Nextclade dataset filenames expected by viralQC."""

from argparse import ArgumentParser
import json
from pathlib import Path


def normalize_dataset_dir(out_dir: Path) -> None:
    """
    Link or stub files so a downloaded dataset matches viralQC layout.

    Some community datasets use custom names in pathogen.json (e.g.
    reference_M.fasta instead of reference.fasta). Snakemake rules expect
    the standard names used across the pipeline.
    """
    pathogen = out_dir / "pathogen.json"
    if pathogen.exists():
        files = json.loads(pathogen.read_text()).get("files", {})
        for std, custom in (
            ("reference.fasta", files.get("reference")),
            ("sequences.fasta", files.get("examples")),
            ("tree.json", files.get("treeJson")),
        ):
            if not custom:
                continue
            src = out_dir / custom
            dst = out_dir / std
            if src.exists() and not dst.exists():
                dst.symlink_to(src.name)

    for pattern, std in (
        ("reference*.fasta", "reference.fasta"),
        ("sample*.fasta", "sequences.fasta"),
    ):
        if not (out_dir / std).exists():
            matches = sorted(out_dir.glob(pattern))
            if matches:
                (out_dir / std).symlink_to(matches[0].name)

    if not (out_dir / "tree.json").exists():
        tree_candidates = sorted(out_dir.glob("*.json"))
        tree_candidates = [p for p in tree_candidates if p.name != "pathogen.json"]
        if tree_candidates:
            (out_dir / "tree.json").symlink_to(tree_candidates[0].name)

    for name in (
        "sequences.fasta",
        "tree.json",
        "README.md",
        "CHANGELOG.md",
        "genome_annotation.gff3",
        "pathogen.json",
        "reference.fasta",
    ):
        path = out_dir / name
        if not path.exists():
            path.touch()


def main() -> None:
    parser = ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory where nextclade dataset get wrote files",
    )
    args = parser.parse_args()
    normalize_dataset_dir(args.output_dir)


if __name__ == "__main__":
    main()
