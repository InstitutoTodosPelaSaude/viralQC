import argparse
import logging
import re
import sys
from pathlib import Path


def sanitize_name(name: str) -> str:
    """
    Replace characters that are unsafe in file names with underscores.

    Args:
        name: Original sample header string.

    Returns:
        String with spaces, slashes and special characters replaced by '_'.
    """
    return re.sub(r"[^A-Za-z0-9_.\-]", "_", name)


def load_id_mapping(mapping_file: Path) -> dict:
    """
    Load the numeric ID to original header mapping from a TSV file.

    Args:
        mapping_file: Path to the TSV file with columns: id, original_header.

    Returns:
        Dictionary mapping str(id) to original_header.
    """
    mapping: dict = {}
    with open(mapping_file, "r", encoding="utf-8") as fh:
        for lineno, line in enumerate(fh, start=1):
            line = line.rstrip("\r\n")
            if lineno == 1 and line.lower().startswith("id\t"):
                continue
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                logging.warning(
                    "Skipping malformed line %d in mapping file: %r", lineno, line
                )
                continue
            mapping[parts[0].strip()] = parts[1].strip()
    return mapping


def parse_tbl_blocks(tbl_file: Path):
    """
    Iterate over per-sample blocks in a nextclade multi-sample TBL file.

    Each block starts with a '>Feature <seqid>' header line and contains
    all annotation lines for that sample.

    Args:
        tbl_file: Path to the nextclade TBL file.

    Yields:
        Tuple of (sample_id, lines) where sample_id is the sequence
        identifier string and lines is the list of raw text lines for
        that sample's block (including the '>Feature' line).
    """
    current_id: str | None = None
    current_lines: list = []

    with open(tbl_file, "r", encoding="utf-8") as fh:
        for raw_line in fh:
            line = raw_line.rstrip("\r\n")

            if line.startswith(">Feature "):
                if current_id is not None:
                    yield (current_id, current_lines)
                current_id = line.split(" ", 1)[1].strip()
                current_lines = [line]
            elif current_id is not None:
                current_lines.append(line)

        if current_id is not None and current_lines:
            yield (current_id, current_lines)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Split a nextclade multi-sample TBL file into one TBL file per sample."
    )
    parser.add_argument(
        "--tbl",
        type=Path,
        required=True,
        help="Path to the nextclade TBL file to split.",
    )
    parser.add_argument(
        "--id-mapping",
        type=Path,
        required=True,
        help="TSV file mapping numeric IDs to original sample headers (columns: id, original_header).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory for per-sample TBL files. Default: <tbl_dir>/per_sample/",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging verbosity level.",
    )

    args = parser.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(levelname)s: %(message)s",
        stream=sys.stderr,
    )

    output_dir = args.output_dir if args.output_dir else args.tbl.parent / "per_sample"
    output_dir.mkdir(parents=True, exist_ok=True)

    id_mapping = load_id_mapping(args.id_mapping)

    n_written = 0
    n_skipped = 0

    for sample_id, block_lines in parse_tbl_blocks(args.tbl):
        if sample_id not in id_mapping:
            logging.warning(
                "No mapping found for sample ID '%s' – skipping block.", sample_id
            )
            n_skipped += 1
            continue

        original_header = id_mapping[sample_id]
        out_filepath = output_dir / f"{sample_id}_{sanitize_name(original_header)}.tbl"

        logging.debug(
            "Writing sample '%s' -> '%s' to %s",
            sample_id,
            original_header,
            out_filepath,
        )

        with open(out_filepath, "w", encoding="utf-8") as out_fh:
            for line in block_lines:
                if line.startswith(">Feature "):
                    out_fh.write(f">Feature {original_header}\n")
                else:
                    out_fh.write(line + "\n")

        n_written += 1

    logging.info(
        "Done. Written: %d sample TBL files. Skipped (no mapping): %d.",
        n_written,
        n_skipped,
    )
