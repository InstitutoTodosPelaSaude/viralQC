import argparse
import logging
import os
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


def collect_global_headers(gff_file: Path) -> list:
    """
    Collect global pragma lines that appear before the first sequence block.

    Args:
        gff_file: Path to the nextclade GFF file.

    Returns:
        List of header lines (e.g. '##gff-version 3') preceding the first
        ##sequence-region directive.
    """
    headers = []
    with open(gff_file, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith("##sequence-region"):
                break
            headers.append(line)
    return headers


def parse_gff_blocks(gff_file: Path):
    """
    Iterate over per-sample blocks in a nextclade multi-sample GFF file.

    Each block starts with a ##sequence-region pragma and contains all
    annotation lines for that sample.

    Args:
        gff_file: Path to the nextclade GFF file.

    Yields:
        Tuple of (sample_id, lines) where sample_id is the sequence
        identifier string and lines is the list of raw text lines for
        that sample's block (including the ##sequence-region line).
    """
    current_id: str | None = None
    current_lines: list = []

    with open(gff_file, "r", encoding="utf-8") as fh:
        for raw_line in fh:
            line = raw_line.rstrip("\r\n")

            if line.startswith("##sequence-region"):
                if current_id is not None:
                    yield (current_id, current_lines)
                    current_lines = []
                parts = line.split()
                current_id = parts[1] if len(parts) >= 2 else "unknown"
                current_lines = [line]
            elif current_id is not None:
                current_lines.append(line)

        if current_id is not None and current_lines:
            yield (current_id, current_lines)


def replace_id_in_line(line: str, old_id: str, new_id: str) -> str:
    """
    Replace the numeric sequence ID with the original sample header in a GFF line.

    Handles ##sequence-region pragmas (column 2) and data lines (column 1
    and bare attribute values such as ID=<id> or Name=<id>).

    Args:
        line: A single line from the GFF file (without trailing newline).
        old_id: Numeric ID to replace.
        new_id: Original sample header to substitute.

    Returns:
        Line with old_id replaced by new_id wherever it appears as a
        sequence identifier.
    """
    if line.startswith("##sequence-region"):
        parts = line.split(" ", 3)
        if len(parts) >= 2 and parts[1] == old_id:
            parts[1] = new_id
        return " ".join(parts)

    if line.startswith("#") or not line:
        return line

    cols = line.split("\t")
    if len(cols) < 9:
        return line

    if cols[0] == old_id:
        cols[0] = new_id

    cols[8] = re.sub(r"(?<==)" + re.escape(old_id) + r"(?=[;,\s]|$)", new_id, cols[8])

    return "\t".join(cols)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Split a nextclade multi-sample GFF file into one GFF file per sample."
    )
    parser.add_argument(
        "--gff",
        type=Path,
        required=True,
        help="Path to the nextclade GFF file to split.",
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
        help="Directory for per-sample GFF files. Default: <gff_dir>/per_sample/",
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

    output_dir = args.output_dir if args.output_dir else args.gff.parent / "per_sample"
    output_dir.mkdir(parents=True, exist_ok=True)

    id_mapping = load_id_mapping(args.id_mapping)
    global_headers = collect_global_headers(args.gff)

    n_written = 0
    n_skipped = 0

    for sample_id, block_lines in parse_gff_blocks(args.gff):
        if sample_id not in id_mapping:
            logging.warning(
                "No mapping found for sample ID '%s' – skipping block.", sample_id
            )
            n_skipped += 1
            continue

        original_header = id_mapping[sample_id]
        out_filepath = output_dir / f"{sample_id}_{sanitize_name(original_header)}.gff"

        logging.debug(
            "Writing sample '%s' -> '%s' to %s",
            sample_id,
            original_header,
            out_filepath,
        )

        with open(out_filepath, "w", encoding="utf-8") as out_fh:
            for hdr in global_headers:
                out_fh.write(hdr + "\n")
            for line in block_lines:
                out_fh.write(
                    replace_id_in_line(line, sample_id, original_header) + "\n"
                )

        n_written += 1

    logging.info(
        "Done. Written: %d sample GFF files. Skipped (no mapping): %d.",
        n_written,
        n_skipped,
    )
