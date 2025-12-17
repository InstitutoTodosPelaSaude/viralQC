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
    subprocess.check_call(cmd)


def main():
    args = parse_args()
    has_spaces = False
    with open(args.input, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            if " " in record.description or " " in record.id:
                pass

    records = list(SeqIO.parse(args.input, "fasta"))
    needs_renaming = False
    for rec in records:
        if rec.description and rec.description != rec.id:
            needs_renaming = True
            break
        if " " in rec.id:
            needs_renaming = True
            break

    if needs_renaming:
        temp_fasta = args.input + ".tmp.fasta"
        mapping_tsv = args.input + ".mapping.tsv"

        mapping = []
        renamed_records = []

        for idx, rec in enumerate(records, 1):
            new_id = str(idx)
            original_name = rec.description
            mapping.append({"id": new_id, "original_name": original_name})

            rec.id = new_id
            rec.description = ""
            renamed_records.append(rec)

        SeqIO.write(renamed_records, temp_fasta, "fasta")
        pd.DataFrame(mapping).to_csv(mapping_tsv, sep="\t", index=False, header=False)

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

        if os.path.exists(args.output) and os.path.getsize(args.output) > 0:
            try:
                df = pd.read_csv(args.output, sep="\t", header=None)
                map_df = pd.read_csv(
                    mapping_tsv,
                    sep="\t",
                    header=None,
                    names=["id", "original_name"],
                    dtype=str,
                )
                id_map = dict(zip(map_df["id"], map_df["original_name"]))
                df[0] = df[0].astype(str).map(id_map)
                df.to_csv(args.output, sep="\t", index=False, header=False)

            except pd.errors.EmptyDataError:
                pass

        if os.path.exists(temp_fasta):
            os.remove(temp_fasta)
        if os.path.exists(mapping_tsv):
            os.remove(mapping_tsv)

    else:
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
