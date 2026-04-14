"""Unit tests for viralqc/scripts/python/*.py — all testable pure functions."""

import csv
import io
import json
import textwrap
from pathlib import Path

import pytest

# ── validate_fasta ────────────────────────────────────────────────────────────

from viralqc.scripts.python.validate_fasta import validate_fasta_file


class TestValidateFasta:
    def _write(self, tmp_path, content, name="test.fasta"):
        p = tmp_path / name
        p.write_text(content, encoding="utf-8")
        return p

    def test_valid_single_sequence(self, tmp_path):
        p = self._write(tmp_path, ">S001\nACGTACGT\n")
        err, total = validate_fasta_file(p)
        assert err is None
        assert total == 1

    def test_valid_multiple_sequences(self, tmp_path):
        content = ">S001\nACGT\n>S002\nTTTT\n>S003\nCCCC\n"
        p = self._write(tmp_path, content)
        err, total = validate_fasta_file(p)
        assert err is None
        assert total == 3

    def test_duplicate_header_detected(self, tmp_path):
        p = self._write(tmp_path, ">S001\nACGT\n>S001\nTTTT\n")
        err, _ = validate_fasta_file(p)
        assert err is not None
        assert "Duplicate" in err

    def test_empty_header_detected(self, tmp_path):
        p = self._write(tmp_path, ">\nACGT\n")
        err, _ = validate_fasta_file(p)
        assert err is not None
        assert "Empty header" in err

    def test_header_without_sequence_detected(self, tmp_path):
        p = self._write(tmp_path, ">S001\n>S002\nACGT\n")
        err, _ = validate_fasta_file(p)
        assert err is not None
        assert "without sequence" in err

    def test_empty_file(self, tmp_path):
        p = self._write(tmp_path, "")
        err, total = validate_fasta_file(p)
        assert err is None
        assert total == 0

    def test_trailing_newlines_ignored(self, tmp_path):
        p = self._write(tmp_path, ">S001\nACGT\n\n\n")
        err, total = validate_fasta_file(p)
        assert err is None
        assert total == 1


# ── reorder_cds ───────────────────────────────────────────────────────────────

from viralqc.scripts.python.reorder_cds import (
    parse_cds_coverage,
    reorder_cds_coverage,
    read_gff_gene_order,
    process_nextclade_tsv,
)


class TestParseCdsCoverage:
    def test_normal(self):
        result = parse_cds_coverage("E:1,NS1:0.9,NS2:0.5")
        assert result == {"E": "1", "NS1": "0.9", "NS2": "0.5"}

    def test_empty_string(self):
        assert parse_cds_coverage("") == {}

    def test_whitespace_only(self):
        assert parse_cds_coverage("   ") == {}

    def test_no_colon(self):
        result = parse_cds_coverage("E1NS1")
        assert result == {}

    def test_single_colon(self):
        result = parse_cds_coverage("E:1")
        assert result == {"E": "1"}


class TestReorderCdsCoverage:
    def test_reorders(self):
        cds_dict = {"NS1": "0.9", "E": "1", "NS2": "0.5"}
        gene_order = ["E", "NS1", "NS2"]
        result = reorder_cds_coverage(cds_dict, gene_order)
        assert result == "E:1,NS1:0.9,NS2:0.5"

    def test_missing_gene_filled_with_zero(self):
        cds_dict = {"E": "1"}
        gene_order = ["E", "NS1"]
        result = reorder_cds_coverage(cds_dict, gene_order)
        assert "NS1:0.0" in result

    def test_empty_order(self):
        result = reorder_cds_coverage({"E": "1"}, [])
        assert result == ""


class TestReadGffGeneOrder:
    def test_extracts_genes_in_order(self, tmp_path):
        gff = tmp_path / "test.gff"
        gff.write_text(
            "##gff-version 3\n"
            "chr1\t.\tgene\t100\t500\t.\t+\t.\tgene_name=E\n"
            "chr1\t.\tgene\t10\t90\t.\t+\t.\tgene_name=NS1\n",
            encoding="utf-8",
        )
        order = read_gff_gene_order(gff)
        assert order == ["NS1", "E"]  # sorted by start position

    def test_skips_comment_lines(self, tmp_path):
        gff = tmp_path / "test.gff"
        gff.write_text("# comment\nchr1\t.\tgene\t1\t100\t.\t+\t.\tgene_name=E\n")
        assert read_gff_gene_order(gff) == ["E"]

    def test_returns_empty_for_no_genes(self, tmp_path):
        gff = tmp_path / "empty.gff"
        gff.write_text("# just a comment\n")
        assert read_gff_gene_order(gff) == []


class TestProcessNextcladeTsv:
    def test_rewrites_cds_coverage(self, tmp_path):
        gff = tmp_path / "test.gff"
        gff.write_text(
            "chr1\t.\tgene\t1\t100\t.\t+\t.\tgene_name=E\n"
            "chr1\t.\tgene\t200\t300\t.\t+\t.\tgene_name=NS1\n"
        )
        tsv_in = tmp_path / "in.tsv"
        tsv_in.write_text("seqName\tcdsCoverage\nS001\tNS1:0.9,E:1.0\n")
        tsv_out = tmp_path / "out.tsv"
        process_nextclade_tsv(tsv_in, gff, tsv_out)
        content = tsv_out.read_text()
        # E comes before NS1 in GFF order
        assert "E:1.0,NS1:0.9" in content


# ── split_tbl_by_sample ────────────────────────────────────────────────────────

from viralqc.scripts.python.split_tbl_by_sample import (
    sanitize_name,
    load_id_mapping,
    parse_tbl_blocks,
)


class TestSanitizeName:
    def test_safe_name_unchanged(self):
        assert sanitize_name("Sample_001") == "Sample_001"

    def test_space_replaced(self):
        assert sanitize_name("sample 001") == "sample_001"

    def test_slash_replaced(self):
        assert sanitize_name("sample/001") == "sample_001"


class TestLoadIdMapping:
    def test_loads_mapping(self, tmp_path):
        mapping_file = tmp_path / "mapping.tsv"
        mapping_file.write_text("id\toriginal_header\n1\tS001\n2\tS002\n")
        result = load_id_mapping(mapping_file)
        assert result["1"] == "S001"
        assert result["2"] == "S002"

    def test_skips_header_row(self, tmp_path):
        p = tmp_path / "mapping.tsv"
        p.write_text("id\toriginal_header\n1\tS001\n")
        result = load_id_mapping(p)
        assert "id" not in result

    def test_handles_no_header(self, tmp_path):
        p = tmp_path / "mapping.tsv"
        p.write_text("1\tS001\n")
        result = load_id_mapping(p)
        assert result["1"] == "S001"


class TestParseTblBlocks:
    def test_yields_blocks(self, tmp_path):
        tbl = tmp_path / "test.tbl"
        tbl.write_text(">Feature 1\n1\t100\tgene\n>Feature 2\n200\t300\tgene\n")
        blocks = list(parse_tbl_blocks(tbl))
        assert len(blocks) == 2
        assert blocks[0][0] == "1"
        assert blocks[1][0] == "2"

    def test_block_includes_feature_line(self, tmp_path):
        tbl = tmp_path / "test.tbl"
        tbl.write_text(">Feature abc\n1\t50\tCDS\n")
        blocks = list(parse_tbl_blocks(tbl))
        assert ">Feature abc" in blocks[0][1][0]

    def test_empty_file_yields_nothing(self, tmp_path):
        tbl = tmp_path / "empty.tbl"
        tbl.write_text("")
        assert list(parse_tbl_blocks(tbl)) == []


# ── format_nextclade_sort ─────────────────────────────────────────────────────

from viralqc.scripts.python.format_nextclade_sort import (
    create_fasta_path,
    map_datasets_to_local_paths,
    format_nextclade_output,
    write_unmapped_sequences,
)


class TestCreateFastaPath:
    def test_normal(self, tmp_path):
        tsv = tmp_path / "results.tsv"
        result = create_fasta_path("dataset_A", tsv, "sequences.fa")
        assert result == tsv.parent / "dataset_A" / "sequences.fa"

    def test_nan_returns_none(self, tmp_path):
        import numpy as np

        tsv = tmp_path / "results.tsv"
        assert create_fasta_path(float("nan"), tsv, "sequences.fa") is None


class TestMapDatasetsToLocalPaths:
    def test_maps_nextclade_datasets(self, tmp_path):
        config = tmp_path / "config.yaml"
        config.write_text(
            "nextclade_data:\n"
            "  sars_cov2:\n"
            "    dataset: nextstrain/sars-cov-2\n"
            "github: {}\n"
        )
        datasets_path = tmp_path / "datasets"
        result = map_datasets_to_local_paths(datasets_path, config)
        assert "nextstrain/sars-cov-2" in result
        assert result["nextstrain/sars-cov-2"] == datasets_path / "sars_cov2"


class TestWriteUnmappedSequences:
    def _make_df(self, rows):
        import pandas as pd

        return pd.DataFrame(rows)

    def test_writes_unmapped(self, tmp_path):
        import pandas as pd

        df = pd.DataFrame(
            [
                {"seqName": "S001", "dataset": "known", "localDataset": tmp_path},
                {"seqName": "S002", "dataset": float("nan"), "localDataset": None},
            ]
        )
        write_unmapped_sequences(df, tmp_path)
        content = (tmp_path / "unmapped_sequences.txt").read_text()
        assert "S002" in content
        assert "S001" not in content

    def test_writes_empty_file_when_all_mapped(self, tmp_path):
        import pandas as pd

        df = pd.DataFrame(
            [
                {"seqName": "S001", "dataset": "ds", "localDataset": tmp_path},
            ]
        )
        write_unmapped_sequences(df, tmp_path)
        content = (tmp_path / "unmapped_sequences.txt").read_text()
        assert content == ""


# ── jsonl_to_gff helpers ──────────────────────────────────────────────────────

from viralqc.scripts.python.jsonl_to_gff import (
    parse_fasta_lengths,
    clean_cds_name,
)


class TestParseFastaLengths:
    def test_single_sequence(self, tmp_path):
        p = tmp_path / "seqs.fasta"
        p.write_text(">acc1\nACGTACGT\n")
        result = parse_fasta_lengths(str(p))
        assert result["acc1"] == 8

    def test_multiple_sequences(self, tmp_path):
        p = tmp_path / "seqs.fasta"
        p.write_text(">acc1\nAAAA\n>acc2\nCCCCCC\n")
        result = parse_fasta_lengths(str(p))
        assert result["acc1"] == 4
        assert result["acc2"] == 6


class TestCleanCdsName:
    def test_spaces_replaced(self):
        result = clean_cds_name("my gene")
        assert " " not in result

    def test_truncated_to_20(self):
        long_name = "A" * 30
        result = clean_cds_name(long_name)
        assert len(result) <= 20

    def test_semicolon_truncation(self):
        result = clean_cds_name("gene_A; isoform_B")
        assert ";" not in result

    def test_colon_truncation(self):
        result = clean_cds_name("gene_A: details")
        assert ":" not in result
