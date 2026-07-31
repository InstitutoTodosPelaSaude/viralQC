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


# ── itemize_sequences ─────────────────────────────────────────────────────────

from viralqc.scripts.python.itemize_sequences import (
    strip_gaps,
    unusable_reason,
    itemize_sequences,
)


class TestStripGaps:
    def test_removes_hyphen_gaps(self):
        assert strip_gaps("ACGT----ACGT") == "ACGTACGT"

    def test_removes_dot_gaps(self):
        assert strip_gaps("ACGT....ACGT") == "ACGTACGT"

    def test_removes_mixed_gap_characters(self):
        assert strip_gaps("AC-GT.ACGT-") == "ACGTACGT"

    def test_preserves_case(self):
        assert strip_gaps("acgt--ACGT") == "acgtACGT"

    def test_leaves_gapless_sequence_untouched(self):
        assert strip_gaps("ACGTNNNN") == "ACGTNNNN"

    def test_all_gaps_becomes_empty(self):
        assert strip_gaps("--------") == ""


class TestUnusableReason:
    def test_normal_sequence_is_usable(self):
        assert unusable_reason("ACGTACGT") is None

    def test_lowercase_sequence_is_usable(self):
        assert unusable_reason("acgtacgt") is None

    def test_partially_gapped_sequence_is_usable(self):
        assert unusable_reason("ACGT----ACGT") is None

    def test_mostly_n_with_one_base_is_usable(self):
        assert unusable_reason("NNNNNNNACNNNNNN") is None

    def test_rna_uracil_is_usable(self):
        assert unusable_reason("ACGUACGU") is None

    def test_all_hyphen_gaps_is_unusable(self):
        reason = unusable_reason("-" * 100)
        assert reason is not None
        assert "gap" in reason

    def test_all_dot_gaps_is_unusable(self):
        reason = unusable_reason("." * 100)
        assert reason is not None
        assert "gap" in reason

    def test_empty_sequence_is_unusable(self):
        reason = unusable_reason("")
        assert reason is not None
        assert "empty" in reason

    def test_all_n_is_unusable(self):
        reason = unusable_reason("N" * 100)
        assert reason is not None
        assert "ambiguous" in reason

    def test_only_ambiguity_codes_is_unusable(self):
        reason = unusable_reason("RYKMSWRYKM")
        assert reason is not None
        assert "ambiguous" in reason


class TestItemizeSequences:
    def _run(self, tmp_path, content):
        src = tmp_path / "in.fasta"
        src.write_text(content, encoding="utf-8")
        fasta = tmp_path / "out.fasta"
        mapping = tmp_path / "map.tsv"
        excluded = tmp_path / "excluded.tsv"
        kept, skipped = itemize_sequences(src, fasta, mapping, excluded)
        return {
            "kept": kept,
            "skipped": skipped,
            "fasta": fasta.read_text(encoding="utf-8"),
            "mapping": list(csv.DictReader(mapping.open(newline=""), delimiter="\t")),
            "excluded": list(csv.DictReader(excluded.open(newline=""), delimiter="\t")),
        }

    def test_renames_to_sequential_ids(self, tmp_path):
        r = self._run(tmp_path, ">first\nACGT\n>second\nTTTT\n")
        assert r["kept"] == 2
        assert r["skipped"] == 0
        assert [m["id"] for m in r["mapping"]] == ["1", "2"]
        assert [m["original_header"] for m in r["mapping"]] == ["first", "second"]
        assert r["excluded"] == []

    def test_gaps_stripped_from_written_fasta(self, tmp_path):
        r = self._run(tmp_path, ">s1\nACGT----ACGT\n>s2\nAC..GT\n")
        assert "-" not in r["fasta"]
        assert "." not in r["fasta"]
        assert "ACGTACGT" in r["fasta"]
        assert "ACGT\n" in r["fasta"] or "ACGT" in r["fasta"]

    def test_all_gap_record_excluded_from_fasta(self, tmp_path):
        r = self._run(tmp_path, ">good\nACGTACGT\n>bad\n" + "-" * 60 + "\n")
        assert r["kept"] == 1
        assert r["skipped"] == 1
        assert ">1" in r["fasta"]
        assert ">2" not in r["fasta"]

    def test_excluded_record_still_present_in_mapping(self, tmp_path):
        r = self._run(tmp_path, ">good\nACGTACGT\n>bad\n" + "-" * 60 + "\n")
        assert [m["original_header"] for m in r["mapping"]] == ["good", "bad"]

    def test_ids_stay_stable_when_a_record_is_excluded(self, tmp_path):
        content = ">a\nACGT\n>b\n" + "-" * 20 + "\n>c\nTTTT\n"
        r = self._run(tmp_path, content)
        ids = {m["original_header"]: m["id"] for m in r["mapping"]}
        assert ids == {"a": "1", "b": "2", "c": "3"}
        assert ">1" in r["fasta"]
        assert ">3" in r["fasta"]
        assert ">2" not in r["fasta"]

    def test_excluded_tsv_records_id_header_reason_and_length(self, tmp_path):
        r = self._run(tmp_path, ">good\nACGT\n>bad\n" + "-" * 60 + "\n")
        assert len(r["excluded"]) == 1
        row = r["excluded"][0]
        assert row["id"] == "2"
        assert row["original_header"] == "bad"
        assert "gap" in row["reason"]
        assert row["length"] == "60"

    def test_all_n_record_excluded(self, tmp_path):
        r = self._run(tmp_path, ">good\nACGT\n>alln\n" + "N" * 40 + "\n")
        assert r["skipped"] == 1
        assert "ambiguous" in r["excluded"][0]["reason"]

    def test_every_sequence_unusable_yields_empty_fasta(self, tmp_path):
        r = self._run(tmp_path, ">a\n" + "-" * 30 + "\n>b\n" + "N" * 30 + "\n")
        assert r["kept"] == 0
        assert r["skipped"] == 2
        assert r["fasta"].strip() == ""
        assert len(r["excluded"]) == 2

    def test_excluded_tsv_has_header_when_nothing_excluded(self, tmp_path):
        r = self._run(tmp_path, ">a\nACGT\n")
        excluded_text = (tmp_path / "excluded.tsv").read_text(encoding="utf-8")
        assert "original_header" in excluded_text
        assert r["excluded"] == []


# ── post_process_nextclade: excluded (unusable) sequences ─────────────────────

from viralqc.scripts.python.post_process_nextclade import (
    TARGET_COLUMNS,
    UNUSABLE_SEQUENCE_STATUS,
    create_excluded_df,
    _write_csv_tsv_output,
)


class TestCreateExcludedDf:
    def _write_excluded(self, tmp_path, rows):
        p = tmp_path / "excluded_sequences.tsv"
        with p.open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(
                fh,
                fieldnames=["id", "original_header", "reason", "length"],
                delimiter="\t",
            )
            w.writeheader()
            for r in rows:
                w.writerow(r)
        return p

    def test_status_column_is_part_of_schema(self):
        assert "inputSequenceStatus" in TARGET_COLUMNS

    def test_status_column_is_last_so_existing_order_is_preserved(self):
        assert list(TARGET_COLUMNS.keys())[-1] == "inputSequenceStatus"

    def test_returns_one_row_per_excluded_sequence(self, tmp_path):
        p = self._write_excluded(
            tmp_path,
            [
                {
                    "id": "2",
                    "original_header": "bad",
                    "reason": "no usable sequence: only gap characters",
                    "length": "60",
                },
                {
                    "id": "5",
                    "original_header": "worse",
                    "reason": "no usable sequence: only ambiguous bases",
                    "length": "40",
                },
            ],
        )
        df = create_excluded_df(p)
        assert len(df) == 2
        assert sorted(df["seqName"].astype(str).tolist()) == ["2", "5"]

    def test_rows_are_marked_as_qc_failures(self, tmp_path):
        p = self._write_excluded(
            tmp_path,
            [
                {
                    "id": "2",
                    "original_header": "bad",
                    "reason": "no usable sequence: only gap characters",
                    "length": "60",
                }
            ],
        )
        df = create_excluded_df(p)
        row = df.iloc[0]
        assert row["inputSequenceStatus"] == UNUSABLE_SEQUENCE_STATUS
        assert row["genomeQuality"] == "D"
        assert row["virus"] == "Unclassified"
        assert row["virus_species"] == "Unclassified"
        assert row["segment"] == "Unsegmented"
        assert float(row["coverage"]) == 0.0

    def test_carries_every_target_column(self, tmp_path):
        p = self._write_excluded(
            tmp_path,
            [
                {
                    "id": "1",
                    "original_header": "bad",
                    "reason": "no usable sequence: empty record",
                    "length": "0",
                }
            ],
        )
        df = create_excluded_df(p)
        for col in TARGET_COLUMNS:
            assert col in df.columns

    def test_seqname_is_the_sanitized_id_for_later_remapping(self, tmp_path):
        p = self._write_excluded(
            tmp_path,
            [
                {
                    "id": "7",
                    "original_header": "orig header",
                    "reason": "no usable sequence: only gap characters",
                    "length": "9",
                }
            ],
        )
        df = create_excluded_df(p)
        assert df.iloc[0]["seqName"] == "7"

    def test_missing_file_yields_empty_dataframe(self, tmp_path):
        df = create_excluded_df(tmp_path / "does_not_exist.tsv")
        assert df.empty

    def test_header_only_file_yields_empty_dataframe(self, tmp_path):
        p = self._write_excluded(tmp_path, [])
        df = create_excluded_df(p)
        assert df.empty


class TestWriteOutputWithOnlyExcludedSequences:
    """An input whose every sequence is unusable must still produce a valid results file."""

    def _excluded_df(self, tmp_path, ids):
        p = tmp_path / "excluded_sequences.tsv"
        with p.open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(
                fh,
                fieldnames=["id", "original_header", "reason", "length"],
                delimiter="\t",
            )
            w.writeheader()
            for i in ids:
                w.writerow(
                    {
                        "id": str(i),
                        "original_header": f"seq{i}",
                        "reason": "no usable sequence: only gap characters",
                        "length": "60",
                    }
                )
        return create_excluded_df(p)

    def test_tsv_written_with_header_when_there_are_no_other_rows(self, tmp_path):
        out = tmp_path / "results.tsv"
        excluded = self._excluded_df(tmp_path, [1, 2])
        _write_csv_tsv_output(
            iter([]),
            None,
            set(),
            out,
            "tsv",
            id_map={"1": "seqA", "2": "seqB"},
            excluded_df=excluded,
        )
        assert out.exists()
        lines = out.read_text(encoding="utf-8").strip().split("\n")
        header = lines[0].split("\t")
        assert header == list(TARGET_COLUMNS.keys())
        assert len(lines) == 3

    def test_original_headers_restored_in_output(self, tmp_path):
        out = tmp_path / "results.tsv"
        excluded = self._excluded_df(tmp_path, [1])
        _write_csv_tsv_output(
            iter([]),
            None,
            set(),
            out,
            "tsv",
            id_map={"1": "ITpS_real_header"},
            excluded_df=excluded,
        )
        rows = list(csv.DictReader(out.open(newline=""), delimiter="\t"))
        assert rows[0]["seqName"] == "ITpS_real_header"
        assert rows[0]["inputSequenceStatus"] == UNUSABLE_SEQUENCE_STATUS


# ── format_nextclade_sort: empty nextclade sort output ────────────────────────

from viralqc.scripts.python.format_nextclade_sort import (
    read_nextclade_sort_tsv,
    format_nextclade_output,
)


class TestReadNextcladeSortTsv:
    """nextclade sort writes a header-less, zero-byte table when handed a FASTA
    with no records, which happens when every input sequence was excluded for
    holding no usable base."""

    def test_reads_a_normal_table(self, tmp_path):
        p = tmp_path / "viruses.tsv"
        p.write_text("index\tseqName\tdataset\tscore\tnumHits\n0\t1\tdenv2\t0.4\t283\n")
        df = read_nextclade_sort_tsv(p)
        assert len(df) == 1
        assert df.iloc[0]["seqName"] == "1"

    def test_empty_file_yields_empty_frame_with_schema(self, tmp_path):
        p = tmp_path / "viruses.tsv"
        p.write_text("")
        df = read_nextclade_sort_tsv(p)
        assert df.empty
        for col in ("index", "seqName", "dataset", "score", "numHits"):
            assert col in df.columns

    def test_header_only_file_yields_empty_frame(self, tmp_path):
        p = tmp_path / "viruses.tsv"
        p.write_text("index\tseqName\tdataset\tscore\tnumHits\n")
        df = read_nextclade_sort_tsv(p)
        assert df.empty

    def test_format_nextclade_output_survives_two_empty_tables(self, tmp_path):
        a = tmp_path / "viruses.tsv"
        b = tmp_path / "viruses.external_datasets.tsv"
        a.write_text("")
        b.write_text("")
        df = format_nextclade_output(a, b, {})
        assert df.empty
        assert "localDataset" in df.columns
        assert "splittedFasta" in df.columns
