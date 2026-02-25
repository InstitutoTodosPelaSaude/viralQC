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

    try:
        with open(args.output_fasta, "w") as output_handle, open(
            args.output_mapping, "w", newline=""
        ) as map_handle:

            writer = csv.DictWriter(
                map_handle, fieldnames=["id", "original_header"], delimiter="\t"
            )
            writer.writeheader()

            with open(args.input, "r") as input_handle:
                for i, record in enumerate(SeqIO.parse(input_handle, "fasta"), 1):
                    original_id = record.description
                    new_id = str(i)

                    writer.writerow({"id": new_id, "original_header": original_id})

                    record.id = new_id
                    record.description = ""
                    SeqIO.write(record, output_handle, "fasta")

    except Exception as e:
        print(f"Error processing sequences: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
