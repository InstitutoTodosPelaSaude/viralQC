"""Unit tests for viralqc.core.ncbi_submission."""

import json
import textwrap
from pathlib import Path

import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from viralqc.core.ncbi_submission import (
    NCBI_BATCH_SIZE,
    _batch_path,
    _n_fraction,
    _serotype_fields,
    copy_tbl,
    deduplicate_case_insensitive,
    filter_sequences,
    is_plain_header_virus,
    load_fasta,
    load_metadata,
    load_results,
    sanitize_seq_name,
    write_dropped_table,
    write_fasta,
    write_rename_log,
    write_submission_metadata,
)
from tests.conftest import GOOD_SEQ, BAD_SEQ, make_fasta, make_metadata_csv, make_results_tsv


# ── sanitize_seq_name ────────────────────────────────────────────────────────


class TestSanitizeSeqName:
    def test_clean_name_unchanged(self):
        name = "sample_001"
        out, changed = sanitize_seq_name(name)
        assert out == name
        assert not changed

    def test_pipe_replaced(self):
        out, changed = sanitize_seq_name("sample|001")
        assert "|" not in out
        assert changed

    def test_non_ascii_replaced(self):
        out, changed = sanitize_seq_name("sampleé001")
        assert "é" not in out
        assert changed

    def test_organism_tag_preserved(self):
        name = "sample/1 [Organism=Dengue virus]"
        out, changed = sanitize_seq_name(name)
        assert "[Organism=Dengue virus]" in out
        assert "/" not in out.split("[")[0]

    def test_whitespace_replaced(self):
        out, changed = sanitize_seq_name("sample 001")
        assert " " not in out
        assert changed


# ── _n_fraction ─────────────────────────────────────────────────────────────


class TestNFraction:
    def test_all_n(self):
        assert _n_fraction("NNNN") == 1.0

    def test_no_n(self):
        assert _n_fraction("ACGT") == 0.0

    def test_half_n(self):
        assert _n_fraction("NNACGT") == pytest.approx(2 / 6)

    def test_empty(self):
        assert _n_fraction("") == 1.0


# ── load_results ─────────────────────────────────────────────────────────────


class TestLoadResults:
    def test_load_tsv(self, tmp_path):
        p = make_results_tsv(tmp_path, [{"seqName": "S001", "virus": "SARS-CoV-2"}])
        df = load_results(p)
        assert "seqName" in df.columns
        assert df["seqName"].iloc[0] == "S001"

    def test_load_csv(self, tmp_path):
        p = tmp_path / "results.csv"
        df = pd.DataFrame([{"seqName": "S001", "virus": "SARS-CoV-2"}])
        df.to_csv(p, sep=";", index=False)
        result = load_results(p)
        assert result["seqName"].iloc[0] == "S001"

    def test_load_json(self, tmp_path):
        p = tmp_path / "results.json"
        df = pd.DataFrame([{"seqName": "S001", "virus": "Dengue virus type 1"}])
        df.to_json(p, orient="table", index=False)
        result = load_results(p)
        assert result["seqName"].iloc[0] == "S001"

    def test_unsupported_extension_raises(self, tmp_path):
        p = tmp_path / "results.xlsx"
        p.write_text("dummy")
        with pytest.raises(ValueError, match="Unsupported"):
            load_results(p)


# ── load_fasta ───────────────────────────────────────────────────────────────


class TestLoadFasta:
    def test_loads_records(self, tmp_path):
        p = make_fasta(tmp_path, "seqs.fasta", [("S001", GOOD_SEQ), ("S002", GOOD_SEQ)])
        records = load_fasta(p)
        assert "S001" in records
        assert "S002" in records
        assert str(records["S001"].seq) == GOOD_SEQ

    def test_empty_fasta(self, tmp_path):
        p = tmp_path / "empty.fasta"
        p.write_text("")
        records = load_fasta(p)
        assert records == {}


# ── filter_sequences ─────────────────────────────────────────────────────────


class TestFilterSequences:
    def _make_df(self, seq_names):
        return pd.DataFrame({"seqName": seq_names})

    def test_passing_sequences_kept(self, tmp_path):
        fasta = {"S001": SeqRecord(Seq(GOOD_SEQ), id="S001")}
        df = self._make_df(["S001"])
        kept, dropped = filter_sequences(fasta, df)
        assert len(kept) == 1
        assert len(dropped) == 0

    def test_high_n_content_dropped(self):
        fasta = {"S001": SeqRecord(Seq(BAD_SEQ), id="S001")}
        df = self._make_df(["S001"])
        kept, dropped = filter_sequences(fasta, df)
        assert len(kept) == 0
        assert dropped[0]["seqName"] == "S001"
        assert "N content" in dropped[0]["reason"]

    def test_too_short_dropped(self):
        fasta = {"S001": SeqRecord(Seq("ACGT"), id="S001")}
        df = self._make_df(["S001"])
        kept, dropped = filter_sequences(fasta, df, min_length=150)
        assert len(kept) == 0
        assert "shorter" in dropped[0]["reason"]

    def test_mixed_sequences(self):
        fasta = {
            "good": SeqRecord(Seq(GOOD_SEQ), id="good"),
            "bad": SeqRecord(Seq(BAD_SEQ), id="bad"),
        }
        df = self._make_df(["good", "bad"])
        kept, dropped = filter_sequences(fasta, df)
        assert len(kept) == 1
        assert kept[0].id == "good"
        assert len(dropped) == 1

    def test_custom_thresholds(self):
        # 30% N (below 50% default but above custom 25%)
        seq_30n = "N" * 30 + "A" * 70
        fasta = {"S001": SeqRecord(Seq(seq_30n), id="S001")}
        df = self._make_df(["S001"])
        kept, dropped = filter_sequences(fasta, df, max_n_fraction=0.25)
        assert len(kept) == 0


# ── write_fasta ──────────────────────────────────────────────────────────────


class TestWriteFasta:
    def _records(self, n=3):
        return [SeqRecord(Seq(GOOD_SEQ), id=f"S{i:03d}") for i in range(1, n + 1)]

    def test_single_file_written(self, tmp_path):
        records = self._records(2)
        out = tmp_path / "out.fasta"
        write_fasta(records, out)
        assert out.exists()
        content = out.read_text()
        assert ">S001" in content
        assert ">S002" in content

    def test_custom_header_fn(self, tmp_path):
        records = self._records(1)
        out = tmp_path / "out.fasta"
        write_fasta(records, out, header_fn=lambda r: f"{r.id}_custom")
        assert "S001_custom" in out.read_text()

    def test_batch_split_above_limit(self, tmp_path):
        n = NCBI_BATCH_SIZE + 5
        records = [SeqRecord(Seq("A" * 200), id=f"S{i:05d}") for i in range(n)]
        out = tmp_path / "sequences.fasta"
        write_fasta(records, out)
        assert not out.exists(), "single file should not exist when batched"
        assert (tmp_path / "sequences.1.fasta").exists()
        assert (tmp_path / "sequences.2.fasta").exists()

    def test_rename_log_populated(self, tmp_path):
        records = [SeqRecord(Seq(GOOD_SEQ), id="sample|with|pipes")]
        out = tmp_path / "out.fasta"
        rename_log = []
        write_fasta(records, out, rename_log=rename_log)
        assert len(rename_log) == 1
        assert "sample|with|pipes" in rename_log[0]["original"]

    def test_clean_name_not_logged(self, tmp_path):
        records = [SeqRecord(Seq(GOOD_SEQ), id="clean_name")]
        out = tmp_path / "out.fasta"
        rename_log = []
        write_fasta(records, out, rename_log=rename_log)
        assert rename_log == []


# ── _batch_path ──────────────────────────────────────────────────────────────


class TestBatchPath:
    def test_fasta(self):
        p = Path("out/sequences.fasta")
        assert _batch_path(p, 1) == Path("out/sequences.1.fasta")
        assert _batch_path(p, 2) == Path("out/sequences.2.fasta")

    def test_csv(self):
        p = Path("out/metadata.csv")
        assert _batch_path(p, 1) == Path("out/metadata.1.csv")


# ── write_dropped_table / write_rename_log ───────────────────────────────────


class TestWriteHelpers:
    def test_dropped_table_written(self, tmp_path):
        dropped = [{"seqName": "S001", "reason": "too short", "length": 10, "n_pct": 0}]
        out = tmp_path / "dropped.tsv"
        write_dropped_table(dropped, out)
        assert out.exists()
        df = pd.read_csv(out, sep="\t")
        assert df.iloc[0]["seqName"] == "S001"

    def test_dropped_table_empty_deletes_existing(self, tmp_path):
        out = tmp_path / "dropped.tsv"
        out.write_text("dummy")
        write_dropped_table([], out)
        assert not out.exists()

    def test_rename_log_written(self, tmp_path):
        log = [{"original": "a|b", "sanitized": "a_b"}]
        out = tmp_path / "rename.tsv"
        write_rename_log(log, out)
        assert out.exists()

    def test_rename_log_empty_deletes_existing(self, tmp_path):
        out = tmp_path / "rename.tsv"
        out.write_text("dummy")
        write_rename_log([], out)
        assert not out.exists()


# ── copy_tbl ─────────────────────────────────────────────────────────────────


class TestCopyTbl:
    def test_copies_file(self, tmp_path):
        src = tmp_path / "sample.tbl"
        src.write_text("some tbl content")
        dest = tmp_path / "dest"
        result = copy_tbl(str(src), dest)
        assert result is True
        assert (dest / "sample.tbl").exists()

    def test_empty_path_returns_false(self, tmp_path):
        assert copy_tbl("", tmp_path / "dest") is False

    def test_nan_returns_false(self, tmp_path):
        assert copy_tbl(float("nan"), tmp_path / "dest") is False  # type: ignore

    def test_missing_file_returns_false(self, tmp_path):
        assert copy_tbl("/nonexistent/path.tbl", tmp_path / "dest") is False


# ── is_plain_header_virus ────────────────────────────────────────────────────


class TestIsPlainHeaderVirus:
    @pytest.mark.parametrize("virus", [
        "SARS-CoV-2",
        "Dengue virus type 1",
        "Dengue virus type 4",
        "Influenza A H1N1",
        "Influenza A H3N2",
        "Influenza B virus (B/Lee/1940)",
        "Influenza C virus (C/Ann Arbor/1/50)",
        "Norovirus GI",
        "Norovirus GII",
    ])
    def test_standard_viruses(self, virus):
        assert is_plain_header_virus(virus) is True

    @pytest.mark.parametrize("virus", [
        "Oropouche virus",
        "Respiratory syncytial virus A",
        "Rhinovirus A",
    ])
    def test_custom_viruses(self, virus):
        assert is_plain_header_virus(virus) is False

    def test_influenza_by_regex(self):
        # Even unlisted subtypes should be plain-header viruses
        assert is_plain_header_virus("Influenza A H5N1") is True


# ── _serotype_fields ─────────────────────────────────────────────────────────


class TestSerotypeFields:
    def test_sars_cov2_returns_empty(self):
        assert _serotype_fields("SARS-CoV-2", "21L") == {}

    def test_dengue_extracts_genotype(self):
        fields = _serotype_fields("Dengue virus type 2", "some_clade")
        assert fields["genotype"] == "2"
        assert fields["serotype"] == "some_clade"

    def test_influenza_extracts_serotype(self):
        fields = _serotype_fields("Influenza A H3N2", "")
        assert fields["serotype"] == "H3N2"

    def test_norovirus_returns_genotype(self):
        fields = _serotype_fields("Norovirus GII.17", "")
        assert "genotype" in fields

    def test_unknown_virus_returns_empty_serotype(self):
        fields = _serotype_fields("Unknown virus", "")
        assert fields.get("serotype") == ""


# ── load_metadata ─────────────────────────────────────────────────────────────


class TestLoadMetadata:
    def test_valid_metadata_loads(self, tmp_path):
        rows = [{
            "Sequence_ID": "S001",
            "geo_loc_name": "Brazil",
            "host": "Homo sapiens",
            "isolate": "S001/2024",
            "collection-date": "2024-01-01",
            "isolation-source": "Serum",
        }]
        p = make_metadata_csv(tmp_path, rows)
        df = load_metadata(p, is_standard=True)
        assert len(df) == 1
        assert df["Sequence_ID"].iloc[0] == "S001"

    def test_missing_required_column_raises(self, tmp_path):
        p = make_metadata_csv(tmp_path, [{"Sequence_ID": "S001"}])
        with pytest.raises(ValueError, match="missing required columns"):
            load_metadata(p, is_standard=True)

    def test_long_sequence_id_raises(self, tmp_path):
        rows = [{
            "Sequence_ID": "X" * 25,  # ≥25 chars → invalid
            "geo_loc_name": "Brazil",
            "host": "Homo sapiens",
            "isolate": "iso",
            "collection-date": "2024-01-01",
            "isolation-source": "Serum",
        }]
        p = make_metadata_csv(tmp_path, rows)
        with pytest.raises(ValueError, match="less than 25"):
            load_metadata(p, is_standard=True)

    def test_custom_virus_only_requires_sequence_id(self, tmp_path):
        p = make_metadata_csv(tmp_path, [{"Sequence_ID": "S001"}])
        df = load_metadata(p, is_standard=False)
        assert len(df) == 1


# ── write_submission_metadata ─────────────────────────────────────────────────


class TestWriteSubmissionMetadata:
    def _make_records(self, ids):
        return [SeqRecord(Seq(GOOD_SEQ), id=i) for i in ids]

    def test_standard_csv_written(self, tmp_path, tmp_results, tmp_metadata):
        df = load_results(tmp_results)
        meta_df = load_metadata(tmp_metadata)
        records = self._make_records(["S001"])
        out = tmp_path / "metadata.csv"
        write_submission_metadata(records, df, meta_df, out, is_standard=True)
        assert out.exists()
        result = pd.read_csv(out)
        assert "Sequence_ID" in result.columns
        assert result.iloc[0]["Sequence_ID"] == "S001"

    def test_no_records_skips_write(self, tmp_path, tmp_results, tmp_metadata):
        df = load_results(tmp_results)
        meta_df = load_metadata(tmp_metadata)
        out = tmp_path / "meta_no_records.csv"
        write_submission_metadata([], df, meta_df, out, is_standard=True)
        # The function should return early without creating the output file
        assert not out.exists()


    def test_custom_virus_tsv_written(self, tmp_path):
        results = pd.DataFrame([{
            "seqName": "S001",
            "virus": "Oropouche virus",
            "clade": "",
            "virus_species": "Oropouche virus",
        }])
        meta = pd.DataFrame([{
            "Sequence_ID": "S001",
            "geo_loc_name": "Brazil",
            "host": "Homo sapiens",
            "isolate": "iso/2024",
            "collection-date": "2024-01-01",
            "isolation-source": "Serum",
        }])
        records = self._make_records(["S001"])
        out = tmp_path / "metadata.tsv"
        write_submission_metadata(records, results, meta, out, is_standard=False, organism="Oropouche virus")
        assert out.exists()
        result = pd.read_csv(out, sep="\t")
        assert "Sequence_ID" in result.columns

    def test_standard_optional_cols_passed_through(self, tmp_path, tmp_results):
        """Optional INSDC columns present in the input CSV appear in output."""
        meta = pd.DataFrame([{
            "Sequence_ID": "S001",
            "geo_loc_name": "Brazil",
            "host": "Homo sapiens",
            "isolate": "iso/2024",
            "collection-date": "2024-01-01",
            "isolation-source": "Serum",
            "lat_lon": "15.77 S 47.93 W",
            "strain": "B/2024",
            "note": "acute case",
        }])
        df = load_results(tmp_results)
        records = self._make_records(["S001"])
        out = tmp_path / "metadata.csv"
        write_submission_metadata(records, df, meta, out, is_standard=True)
        result = pd.read_csv(out)
        assert "lat_lon" in result.columns
        assert "strain" in result.columns
        assert "note" in result.columns
        assert result.iloc[0]["lat_lon"] == "15.77 S 47.93 W"

    def test_custom_optional_cols_renamed_and_passed_through(self, tmp_path):
        """Optional columns for custom viruses are renamed to BankIt Title_Case."""
        results = pd.DataFrame([{
            "seqName": "S001",
            "virus": "Oropouche virus",
            "clade": "",
            "virus_species": "Oropouche virus",
        }])
        meta = pd.DataFrame([{
            "Sequence_ID": "S001",
            "geo_loc_name": "Brazil",
            "host": "Homo sapiens",
            "isolate": "iso/2024",
            "collection-date": "2024-01-01",
            "isolation-source": "Serum",
            "lat_lon": "15.77 S 47.93 W",
            "strain": "OROV/2024",
            "note": "febrile",
        }])
        records = self._make_records(["S001"])
        out = tmp_path / "metadata.tsv"
        write_submission_metadata(records, results, meta, out, is_standard=False)
        result = pd.read_csv(out, sep="\t")
        assert "Lat_Lon" in result.columns
        assert "Strain" in result.columns
        assert "Note" in result.columns
        assert result.iloc[0]["Strain"] == "OROV/2024"
        # required columns also present
        assert "Country (geo_loc_name)" in result.columns

    def test_unknown_cols_not_passed_through(self, tmp_path, tmp_results):
        """Columns not in either NCBI modifier list are silently excluded."""
        meta = pd.DataFrame([{
            "Sequence_ID": "S001",
            "geo_loc_name": "Brazil",
            "host": "Homo sapiens",
            "isolate": "iso/2024",
            "collection-date": "2024-01-01",
            "isolation-source": "Serum",
            "foo_bar": "should_be_ignored",
        }])
        df = load_results(tmp_results)
        records = self._make_records(["S001"])
        out = tmp_path / "metadata.csv"
        write_submission_metadata(records, df, meta, out, is_standard=True)
        result = pd.read_csv(out)
        assert "foo_bar" not in result.columns


# ── deduplicate_case_insensitive ──────────────────────────────────────────────


class TestDeduplicateCaseInsensitive:
    def _rec(self, id_):
        return SeqRecord(Seq(GOOD_SEQ), id=id_)

    def test_no_duplicates_unchanged(self):
        records = [self._rec("S001"), self._rec("S002"), self._rec("S003")]
        result = deduplicate_case_insensitive(records)
        assert len(result) == 3

    def test_exact_case_duplicate_kept_first(self):
        records = [self._rec("S001"), self._rec("S001")]
        result = deduplicate_case_insensitive(records)
        assert len(result) == 1
        assert result[0].id == "S001"

    def test_case_insensitive_duplicate_dropped(self):
        # ITpS_DLC_000zN_001 and ITpS_DLc_000zn_001 — same ID, different case
        records = [self._rec("ITpS_DLC_000zN_001"), self._rec("ITpS_DLc_000zn_001")]
        result = deduplicate_case_insensitive(records)
        assert len(result) == 1
        assert result[0].id == "ITpS_DLC_000zN_001"  # first one kept

    def test_dropped_appended_to_log(self):
        records = [self._rec("SEQ1"), self._rec("seq1")]
        log = []
        deduplicate_case_insensitive(records, dropped_log=log)
        assert len(log) == 1
        assert log[0]["seqName"] == "seq1"
        assert "Case-insensitive duplicate" in log[0]["reason"]
        assert "SEQ1" in log[0]["reason"]

    def test_multiple_duplicates_all_logged(self):
        records = [self._rec("A"), self._rec("a"), self._rec("A")]
        log = []
        result = deduplicate_case_insensitive(records, dropped_log=log)
        assert len(result) == 1
        assert len(log) == 2

    def test_custom_header_fn_used(self):
        # header_fn maps rec.id → "PREFIX_" + rec.id; collision on lowercase
        records = [self._rec("s1"), self._rec("S1")]
        header_fn = lambda r: r.id.lower()  # both map to "s1"
        log = []
        result = deduplicate_case_insensitive(records, header_fn=header_fn, dropped_log=log)
        assert len(result) == 1
        assert len(log) == 1

    def test_write_fasta_drops_case_duplicates(self, tmp_path):
        records = [self._rec("ITpS_DLC_000zN_001"), self._rec("ITpS_DLc_000zn_001")]
        out = tmp_path / "out.fasta"
        log = []
        write_fasta(records, out, dropped_log=log)
        content = out.read_text()
        # Only the first ID should appear
        assert "ITpS_DLC_000zN_001" in content
        assert "ITpS_DLc_000zn_001" not in content
        assert len(log) == 1
        assert "Case-insensitive" in log[0]["reason"]
