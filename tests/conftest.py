"""
Shared pytest fixtures and helpers for viralqc test suite.
"""

import io
import textwrap
from pathlib import Path

import pandas as pd
import pytest


# ---------------------------------------------------------------------------
# Generic helpers
# ---------------------------------------------------------------------------

def make_fasta(tmp_path: Path, filename: str, records: list[tuple[str, str]]) -> Path:
    """Write a minimal FASTA file. records is a list of (header, sequence)."""
    p = tmp_path / filename
    lines = []
    for hdr, seq in records:
        lines.append(f">{hdr}")
        lines.append(seq)
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return p


def make_results_tsv(tmp_path: Path, rows: list[dict], filename: str = "results.tsv") -> Path:
    """Write a minimal ViralQC results TSV file."""
    p = tmp_path / filename
    df = pd.DataFrame(rows)
    df.to_csv(p, sep="\t", index=False)
    return p


def make_metadata_csv(tmp_path: Path, rows: list[dict], filename: str = "metadata.csv") -> Path:
    """Write a minimal metadata CSV file."""
    p = tmp_path / filename
    df = pd.DataFrame(rows)
    df.to_csv(p, index=False)
    return p


# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

STANDARD_METADATA_ROW = {
    "Sequence_ID": "S001",
    "geo_loc_name": "Brazil",
    "host": "Homo sapiens",
    "isolate": "S001/2024",
    "collection-date": "2024-01-01",
    "isolation-source": "Serum",
}

GOOD_SEQ = "ACGTACGTACGTACGT" * 20  # 320 nt, 0% N → passes filters
BAD_SEQ = "N" * 300  # 100% N → fails filter


@pytest.fixture()
def tmp_fasta(tmp_path):
    """Return a FASTA file with two passing sequences."""
    return make_fasta(tmp_path, "seqs.fasta", [
        ("S001", GOOD_SEQ),
        ("S002", GOOD_SEQ),
    ])


@pytest.fixture()
def tmp_results(tmp_path):
    """Return a minimal TSV results file for SARS-CoV-2."""
    rows = [
        {"seqName": "S001", "virus": "SARS-CoV-2", "clade": "21L", "virus_species": "SARS-CoV-2"},
        {"seqName": "S002", "virus": "SARS-CoV-2", "clade": "21L", "virus_species": "SARS-CoV-2"},
    ]
    return make_results_tsv(tmp_path, rows)


@pytest.fixture()
def tmp_metadata(tmp_path):
    """Return a minimal metadata CSV with two standard-virus rows."""
    rows = [
        {**STANDARD_METADATA_ROW, "Sequence_ID": "S001"},
        {**STANDARD_METADATA_ROW, "Sequence_ID": "S002"},
    ]
    return make_metadata_csv(tmp_path, rows)
