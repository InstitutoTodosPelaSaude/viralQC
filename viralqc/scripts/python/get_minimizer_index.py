#!/usr/bin/env python3

## This a short version of minimizer.py from nextclade project (https://github.com/nextstrain/nextclade_data/blob/master/scripts/minimizer)
## The import paths have been adjusted to work within viralqc project structure.
## The fasta_read was modified to include dataset information in the SeqRecord annotations to
## allow proper reference naming when building the local minimizer index.

# -- lib.fasta
from typing import List, Union, Iterable
from Bio import SeqIO


def is_iterable(obj):
    return issubclass(type(obj), Iterable)


def fasta_read(filepaths: Union[str, List[str]]):
    if isinstance(filepaths, str) or (not is_iterable(filepaths)):
        filepaths = [filepaths]

    records = []
    for filepath in filepaths:
        for record in SeqIO.parse(filepath, "fasta"):
            record.annotations["dataset"] = filepath.split("/")[-2]
            records.append(record)
    return records


# -- lib.fs
import json
from os.path import dirname
from os import makedirs


def ensure_dir(file_path):
    dir_path = dirname(file_path)
    if len(dir_path.strip()) == 0:
        return
    makedirs(dir_path, exist_ok=True)


def file_write(data, filepath):
    ensure_dir(filepath)
    with open(filepath, "w") as f:
        f.write(f"{data.strip()}\n")


def json_write(obj, filepath, no_sort_keys=False):
    content = json.dumps(obj, indent=2, sort_keys=not no_sort_keys, ensure_ascii=False)
    content += "\n"
    file_write(content, filepath)


# -- lib.minimizer
import numpy as np
import copy
from Bio.SeqRecord import SeqRecord

CUTOFF = 1 << 28
MINIMIZER_ALGO_VERSION = "1"
MINIMIZER_JSON_SCHEMA_VERSION = "3.0.0"
JSON_SCHEMA_URL_MINIMIZER_JSON = "https://raw.githubusercontent.com/nextstrain/nextclade/refs/heads/release/packages/nextclade-schemas/internal-minimizer-index-json.schema.json"


def preprocess_seq(seq: SeqRecord) -> str:
    return str(seq.seq).upper().replace("-", "")


def invertible_hash(x):
    m = (1 << 32) - 1
    x = (~x + (x << 21)) & m
    x = x ^ (x >> 24)
    x = (x + (x << 3) + (x << 8)) & m
    x = x ^ (x >> 14)
    x = (x + (x << 2) + (x << 4)) & m
    x = x ^ (x >> 28)
    x = (x + (x << 31)) & m
    return x


def get_hash(kmer):
    x = 0
    j = 0
    for i, nuc in enumerate(kmer):
        if i % 3 == 2:
            continue
        if nuc not in "ACGT":
            return CUTOFF + 1
        else:
            if nuc in "AC":
                x += 1 << j
            if nuc in "AT":
                x += 1 << (j + 1)
        j += 2

    return invertible_hash(x)


def get_ref_search_minimizers(seq: SeqRecord, k=17):
    seq_str = preprocess_seq(seq)
    minimizers = []
    for i in range(len(seq_str) - k):
        kmer = seq_str[i : i + k]
        mhash = get_hash(kmer)
        if mhash < CUTOFF:
            minimizers.append(mhash)
    return np.unique(minimizers)


def make_ref_search_index(refs):
    minimizers_by_reference = list()
    for name, ref in sorted(refs.items()):
        minimizers = get_ref_search_minimizers(ref)
        minimizers_by_reference.append(
            {
                "minimizers": minimizers,
                "meta": {
                    "name": name,
                    "length": len(ref.seq),
                    "nMinimizers": len(minimizers),
                },
            }
        )

    index = {"minimizers": {}, "references": []}
    for ri, minimizer_set in enumerate(minimizers_by_reference):
        for m in minimizer_set["minimizers"]:
            if m not in index["minimizers"]:
                index["minimizers"][m] = []
            index["minimizers"][m].append(ri)

        index["references"].append(minimizer_set["meta"])

    normalization = np.array(
        [x["length"] / x["nMinimizers"] for x in index["references"]]
    )

    return {
        "$schema": JSON_SCHEMA_URL_MINIMIZER_JSON,
        "schemaVersion": MINIMIZER_JSON_SCHEMA_VERSION,
        "version": MINIMIZER_ALGO_VERSION,
        "params": {
            "k": 17,
            "cutoff": CUTOFF,
        },
        **index,
        "normalization": normalization,
    }


def serialize_ref_search_index(index):
    index = copy.deepcopy(index)
    index["minimizers"] = {str(k): v for k, v in index["minimizers"].items()}
    index["normalization"] = index["normalization"].tolist()
    return index


import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "--input-refs",
        required=True,
        nargs="+",
        help="One or more fasta files with reference sequences",
    )
    parser.add_argument(
        "--output-json", required=False, help="Where to output minimizer index file"
    )

    args = parser.parse_args()

    return args


def main():
    args = parse_args()

    refs = fasta_read(args.input_refs)
    refs = {ref.annotations.get("dataset"): ref for ref in refs}
    index = make_ref_search_index(refs)

    if args.output_json:
        json_write(
            serialize_ref_search_index(index), args.output_json, no_sort_keys=True
        )


if __name__ == "__main__":
    main()
