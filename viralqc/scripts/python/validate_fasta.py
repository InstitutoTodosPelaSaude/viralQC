import argparse
import sys
from pathlib import Path


def validate_fasta_file(fasta_path: Path) -> tuple[str | None, int]:
    """
    Validate a FASTA file. Stops at the first error.

    Args:
        fasta_path: Path to the FASTA file to validate.

    Returns:
        Tuple of (error_message or None, total_sequences)
    """
    seen_headers: set[str] = set()
    current_header: str | None = None
    has_sequence: bool = False
    total_sequences: int = 0

    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if current_header is not None and not has_sequence:
                    return f"Header without sequence: {current_header}", total_sequences

                current_header = line[1:].split()[0] if line[1:].strip() else ""
                if not current_header:
                    return "Empty header found", total_sequences

                if current_header in seen_headers:
                    return f"Duplicate header: {current_header}", total_sequences

                seen_headers.add(current_header)
                total_sequences += 1
                has_sequence = False
            else:
                has_sequence = True

    if current_header is not None and not has_sequence:
        return f"Header without sequence: {current_header}", total_sequences

    return None, total_sequences


def main():
    parser = argparse.ArgumentParser(description="Validate FASTA file.")
    parser.add_argument("--input", type=Path, required=True, help="FASTA file path.")
    parser.add_argument(
        "--output", type=Path, required=True, help="Output report path."
    )
    args = parser.parse_args()

    if not args.input.exists():
        print(f"ERROR: File not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    error, total = validate_fasta_file(args.input)

    with open(args.output, "w") as f:
        if error:
            f.write(f"FAILED: {error}\n")
        else:
            f.write(f"PASSED: {total} sequences validated.\n")

    if error:
        print(f"FASTA validation FAILED: {error}", file=sys.stderr)
        sys.exit(1)
    else:
        print(f"FASTA validation PASSED: {total} sequences.")


if __name__ == "__main__":
    main()
