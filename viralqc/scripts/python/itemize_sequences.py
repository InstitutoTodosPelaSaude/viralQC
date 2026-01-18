#!/usr/bin/env python3
import argparse
import sys
from Bio import SeqIO
import csv


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
    return parser.parse_args()


def main():
    args = parse_args()

    mapping = []
    sequences = []

    try:
        with open(args.input, "r") as handle:
            for i, record in enumerate(SeqIO.parse(handle, "fasta"), 1):
                original_id = record.description
                new_id = str(i)

                mapping.append({"id": new_id, "original_header": original_id})

                record.id = new_id
                record.description = ""
                sequences.append(record)

        if sequences:
            with open(args.output_fasta, "w") as output_handle:
                SeqIO.write(sequences, output_handle, "fasta")
        else:
            open(args.output_fasta, "a").close()

        with open(args.output_mapping, "w", newline="") as tsvfile:
            writer = csv.DictWriter(
                tsvfile, fieldnames=["id", "original_header"], delimiter="\t"
            )
            writer.writeheader()
            for row in mapping:
                writer.writerow(row)

    except Exception as e:
        print(f"Error processing sequences: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
