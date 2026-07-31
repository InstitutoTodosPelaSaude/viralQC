import argparse
import sys
from contextlib import ExitStack
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
import csv

# Alignment gap characters. blastn's CFastaReader rejects the whole query file
# once a line falls below roughly one third plausible residues (measured against a
# real database: ~30% residues aborts, ~34% passes), so both an all-gap record and
# a heavily gapped partial genome abort it. Under `set -euo pipefail` that killed
# the workflow and discarded every other sequence's results. Below that threshold
# the same characters are instead dropped silently by blastn, and handled
# differently again by nextclade, so query coordinates and lengths stop matching
# the input. Strip them here rather than relying on each tool's behaviour.
#
# This is not the only padding convention in the wild ('?' is used for missing
# data and still aborts blastn), which is why the blast rule additionally treats a
# non-zero blastn exit as non-fatal instead of relying on this list being complete.
GAP_CHARACTERS = "-."

# Characters that carry usable information. Ambiguity codes (N, R, Y, ...) are
# deliberately excluded: a record composed only of them cannot be assigned to a
# dataset nor scored, so it is reported as a quality-control failure instead of
# being sent through nextclade and BLAST.
USABLE_BASES = frozenset("ACGTU")

# Reason values written to the `reason` column of excluded_sequences.tsv. Part of
# that file's contract, so they live here rather than inline.
REASON_EMPTY = "no usable sequence: empty record"
REASON_ONLY_GAPS = "no usable sequence: only gap characters"
REASON_ONLY_AMBIGUOUS = "no usable sequence: only ambiguous bases"

_GAP_TRANSLATION = str.maketrans("", "", GAP_CHARACTERS)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Rename sequences in a FASTA file to sequential numeric IDs and output a mapping file."
    )
    parser.add_argument("--input", required=True, help="Input FASTA file")
    parser.add_argument(
        "--output-fasta", required=True, help="Output sanitized FASTA file"
    )
    parser.add_argument(
        "--output-mapping", required=True, help="Output mapping TSV file"
    )
    parser.add_argument(
        "--output-excluded",
        required=True,
        help="Output TSV file listing sequences excluded for having no usable bases",
    )
    parser.add_argument(
        "--output-ungapped",
        required=True,
        help=(
            "Output FASTA with gaps stripped but the ORIGINAL headers kept. Used "
            "downstream wherever coordinates produced by nextclade must be applied "
            "back to the user's sequences."
        ),
    )
    return parser.parse_args()


def strip_gaps(sequence: str) -> str:
    """
    Remove alignment gap characters from a sequence.

    Args:
        sequence: Nucleotide sequence, possibly aligned.

    Returns:
        The sequence with every gap character removed, case preserved.
    """
    return sequence.translate(_GAP_TRANSLATION)


def unusable_reason(sequence: str) -> str | None:
    """
    Determine whether a sequence carries any base that can be analysed.

    Args:
        sequence: Nucleotide sequence as read from the input FASTA.

    Returns:
        A human readable reason when the sequence cannot be analysed, or None
        when it holds at least one usable base.
    """
    if not sequence:
        return REASON_EMPTY

    ungapped = strip_gaps(sequence)
    if not ungapped:
        return REASON_ONLY_GAPS

    if not set(ungapped.upper()) & USABLE_BASES:
        return REASON_ONLY_AMBIGUOUS

    return None


def itemize_sequences(
    input_path: Path,
    output_fasta: Path,
    output_mapping: Path,
    output_excluded: Path,
    output_ungapped: Path | None = None,
) -> tuple[int, int]:
    """
    Rename sequences to sequential numeric IDs, strip alignment gaps and set
    aside records that hold no usable base.

    Args:
        input_path: Input FASTA file.
        output_fasta: Sanitized FASTA file, written without excluded records.
        output_mapping: Mapping TSV file, written with every input record.
        output_excluded: TSV file listing the excluded records and why.
        output_ungapped: Optional FASTA with gaps stripped but original headers
            kept. Every input record is written, so it is a faithful
            gap-stripped copy of the input.

    Returns:
        Tuple of (sequences written, sequences excluded).
    """
    kept = 0
    skipped = 0

    with ExitStack() as stack:
        fasta_handle = stack.enter_context(open(output_fasta, "w"))
        map_handle = stack.enter_context(open(output_mapping, "w", newline=""))
        excluded_handle = stack.enter_context(open(output_excluded, "w", newline=""))
        ungapped_handle = (
            stack.enter_context(open(output_ungapped, "w"))
            if output_ungapped is not None
            else None
        )

        mapping_writer = csv.DictWriter(
            map_handle, fieldnames=["id", "original_header"], delimiter="\t"
        )
        mapping_writer.writeheader()

        excluded_writer = csv.DictWriter(
            excluded_handle,
            fieldnames=["id", "original_header", "reason", "length"],
            delimiter="\t",
        )
        excluded_writer.writeheader()

        input_handle = stack.enter_context(open(input_path, "r"))
        for i, record in enumerate(SeqIO.parse(input_handle, "fasta"), 1):
            original_id = record.description
            new_id = str(i)

            # IDs follow the input order and every record is mapped, even
            # the excluded ones: post-processing needs the mapping to
            # restore the original header of the results row it synthesises
            # for each excluded sequence.
            mapping_writer.writerow({"id": new_id, "original_header": original_id})

            sequence = str(record.seq)
            ungapped = strip_gaps(sequence)

            # Written for every record, under the original header, so that
            # coordinates nextclade derived from the gap-stripped sequence can be
            # applied back to the user's sequences in the same frame.
            if ungapped_handle is not None:
                ungapped_handle.write(f">{original_id}\n{ungapped}\n")

            reason = unusable_reason(sequence)
            if reason is not None:
                excluded_writer.writerow(
                    {
                        "id": new_id,
                        "original_header": original_id,
                        "reason": reason,
                        "length": len(sequence),
                    }
                )
                skipped += 1
                continue

            record.seq = Seq(ungapped)
            record.id = new_id
            record.description = ""
            SeqIO.write(record, fasta_handle, "fasta")
            kept += 1

    return kept, skipped


def main():
    args = parse_args()

    try:
        kept, skipped = itemize_sequences(
            Path(args.input),
            Path(args.output_fasta),
            Path(args.output_mapping),
            Path(args.output_excluded),
            Path(args.output_ungapped),
        )
    except Exception as e:
        print(f"Error processing sequences: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Itemized {kept} sequence(s).")
    if skipped:
        print(
            f"WARNING: {skipped} sequence(s) hold no usable base and were excluded "
            f"from nextclade and BLAST; see {args.output_excluded}. They are still "
            "reported in the results file with "
            "inputSequenceStatus='no usable sequence'.",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
