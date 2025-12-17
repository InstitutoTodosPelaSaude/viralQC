#!/usr/bin/env python3
import sys
import subprocess
import argparse
import pandas as pd
from Bio import SeqIO
import os
import shutil


def parse_args():
    parser = argparse.ArgumentParser(
        description="Wrapper for BLAST to handle sequence headers with spaces."
    )
    parser.add_argument("--input", required=True, help="Input FASTA file")
    parser.add_argument("--db", required=True, help="BLAST database")
    parser.add_argument("--output", required=True, help="Output BLAST results TSV")
    parser.add_argument(
        "--task", required=True, help="BLAST task (e.g., megablast, blastn)"
    )
    parser.add_argument("--evalue", required=True, help="E-value threshold")
    parser.add_argument("--qcov", required=True, help="Query coverage HSP percentage")
    parser.add_argument(
        "--perc_identity", required=True, help="Percent identity threshold"
    )
    parser.add_argument("--threads", required=True, help="Number of threads")
    parser.add_argument(
        "--outfmt", required=True, help="Output format string (cols 1-13)"
    )
    return parser.parse_args()


def run_blast(query, db, output, task, evalue, qcov, perc_identity, threads, outfmt):
    cmd = [
        "blastn",
        "-db",
        db,
        "-query",
        query,
        "-out",
        output,
        "-task",
        task,
        "-evalue",
        evalue,
        "-qcov_hsp_perc",
        qcov,
        "-outfmt",
        outfmt,
        "-max_hsps",
        "1",
        "-max_target_seqs",
        "1",
        "-perc_identity",
        perc_identity,
        "-num_threads",
        threads,
    ]
    # For debugging purposes
    print(f"Running command: {' '.join(cmd)}", file=sys.stderr)
    subprocess.check_call(cmd)


def main():
    args = parse_args()

    # Check for spaces in headers
    has_spaces = False
    with open(args.input, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            if " " in record.description or " " in record.id:
                # BioPython parses just the ID up to the first space if description is used,
                # but we want to check the full header line basically.
                # Actually SeqIO id is just the first part. description is the whole thing.
                # If we rebuild, we want to know if 'id' vs 'description' issues happen.
                # The user says "blast breaks name if it contains space".
                # Usually BLAST uses the first word as ID.
                # We want to preserve the FULL header if possible or at least a unique identifier that maps back.
                pass
            # Just iterating all to be safe? No, let's just do the procedure:
            # Re-read to process.

    # We will ALWAYS rename to be consistent and safe, or check first.
    # The requirement says:
    # "- checar se algum nome no input sequences possui espaço."
    # "- caso exista pelo menos 1 sequencia com espaço no nome, o fluxo deve criar um arquivo de .tsv com 2 colunas."

    records = list(SeqIO.parse(args.input, "fasta"))
    needs_renaming = False
    for rec in records:
        # Check if original header has spaces.
        # SeqIO.parse splits id and description.
        # The full header usually is `>id description`.
        # If there is a description, it implies space was present.
        if rec.description and rec.description != rec.id:
            needs_renaming = True
            break
        if " " in rec.id:  # Should not happen with SeqIO usually but checking
            needs_renaming = True
            break

    if needs_renaming:
        print("Spaces found in headers. Renaming sequences...", file=sys.stderr)

        # Create mapping and temp fasta
        temp_fasta = args.input + ".tmp.fasta"
        mapping_tsv = args.input + ".mapping.tsv"

        mapping = []
        renamed_records = []

        for idx, rec in enumerate(records, 1):
            new_id = str(idx)
            original_name = rec.description  # Use full description as the name
            mapping.append({"id": new_id, "original_name": original_name})

            rec.id = new_id
            rec.description = ""
            renamed_records.append(rec)

        SeqIO.write(renamed_records, temp_fasta, "fasta")

        # Save mapping
        pd.DataFrame(mapping).to_csv(mapping_tsv, sep="\t", index=False, header=False)

        # Run BLAST with temp fasta
        run_blast(
            temp_fasta,
            args.db,
            args.output,
            args.task,
            args.evalue,
            args.qcov,
            args.perc_identity,
            args.threads,
            args.outfmt,
        )

        # Remap results
        # Read the BLAST output
        # If the file is empty (touch was called or no results), handle it
        if os.path.exists(args.output) and os.path.getsize(args.output) > 0:
            # Load blast results. We assume no header likely based on params or specified columns?
            # User specified -outfmt "6 qseqid ..." which has no header by default
            # But we should be careful.
            # The tool call instruction was to output specific columns.
            # Let's read it with the columns specified in Snakemake if possible or just generically.

            # The columns passed in outfmt string are space separated: "6 qseqid qlen ..."
            # So the file is tab separated.
            # Column 1 is qseqid.

            try:
                # We need to read it as CSV/TSV
                df = pd.read_csv(args.output, sep="\t", header=None)

                # Load mapping
                map_df = pd.read_csv(
                    mapping_tsv,
                    sep="\t",
                    header=None,
                    names=["id", "original_name"],
                    dtype=str,
                )
                id_map = dict(zip(map_df["id"], map_df["original_name"]))

                # Replace qseqid (col 0)
                # Ensure col 0 is string to match valid keys
                df[0] = df[0].astype(str).map(id_map)

                # Save back
                df.to_csv(args.output, sep="\t", index=False, header=False)

            except pd.errors.EmptyDataError:
                pass  # Empty file

        # Cleanup
        if os.path.exists(temp_fasta):
            os.remove(temp_fasta)
        if os.path.exists(mapping_tsv):
            # User didn't ask to keep it, but said "create a file". I will keep it or delete?
            # "o fluxo deve criar um arquivo de .tsv com 2 colunas... deve remapear... usando como referencia"
            # It implies it's a temporary intermediate for the remapping.
            # I will delete it to avoid clutter unless requested otherwise.
            os.remove(mapping_tsv)

    else:
        print("No spaces in headers. Running BLAST directly.", file=sys.stderr)
        run_blast(
            args.input,
            args.db,
            args.output,
            args.task,
            args.evalue,
            args.qcov,
            args.perc_identity,
            args.threads,
            args.outfmt,
        )


if __name__ == "__main__":
    main()
