"""Unit tests for viralqc.core.run_analysis.RunAnalysis."""

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from viralqc.core.errors import InvalidOutputFormat
from viralqc.core.run_analysis import RunAnalysis


class TestGetOutputFormat:
    def setup_method(self):
        self.ra = RunAnalysis()

    @pytest.mark.parametrize("ext", ["csv", "tsv", "json"])
    def test_valid_extensions(self, ext):
        fmt = self.ra._get_output_format(f"output.{ext}")
        assert fmt == ext

    def test_invalid_extension_raises(self):
        with pytest.raises(InvalidOutputFormat):
            self.ra._get_output_format("output.xlsx")

    def test_no_extension_raises(self):
        with pytest.raises((InvalidOutputFormat, Exception)):
            self.ra._get_output_format("outputnoext")


class TestRun:
    def setup_method(self):
        self.ra = RunAnalysis()

    def _mock_fasta(self, tmp_path) -> Path:
        p = tmp_path / "seqs.fasta"
        p.write_text(">S001\nACGT\n")
        return p

    @patch("viralqc.core.run_analysis.run_snakemake")
    def test_run_calls_snakemake(self, mock_snakemake, tmp_path):
        """run() must call run_snakemake exactly once with expected keys."""
        fasta = self._mock_fasta(tmp_path)
        mock_snakemake.return_value = MagicMock(success=True)

        self.ra.run(sequences_fasta=fasta, output_file=tmp_path / "results.tsv")

        mock_snakemake.assert_called_once()
        call_kwargs = mock_snakemake.call_args[1]
        assert "config" in call_kwargs
        assert "sequences_fasta" in call_kwargs["config"]

    @patch("viralqc.core.run_analysis.run_snakemake")
    def test_config_contains_expected_keys(self, mock_snakemake, tmp_path):
        fasta = self._mock_fasta(tmp_path)
        mock_snakemake.return_value = MagicMock(success=True)

        self.ra.run(sequences_fasta=fasta, output_file=tmp_path / "results.tsv")

        config = mock_snakemake.call_args[1]["config"]
        for key in (
            "sequences_fasta",
            "output_dir",
            "output_file",
            "output_format",
            "threads",
        ):
            assert key in config, f"Expected key '{key}' in snakemake config"

    @patch("viralqc.core.run_analysis.run_snakemake")
    def test_run_returns_snakemake_response(self, mock_snakemake, tmp_path):
        fasta = self._mock_fasta(tmp_path)
        fake_response = MagicMock(success=True)
        mock_snakemake.return_value = fake_response

        result = self.ra.run(
            sequences_fasta=fasta, output_file=tmp_path / "results.tsv"
        )
        assert result is fake_response

    @patch("viralqc.core.run_analysis.run_snakemake")
    def test_invalid_output_format_raises_before_snakemake(
        self, mock_snakemake, tmp_path
    ):
        fasta = self._mock_fasta(tmp_path)

        with pytest.raises(InvalidOutputFormat):
            self.ra.run(sequences_fasta=fasta, output_file=tmp_path / "results.xlsx")

        mock_snakemake.assert_not_called()

    @patch("viralqc.core.run_analysis.run_snakemake")
    def test_fasta_path_resolved_to_absolute(self, mock_snakemake, tmp_path):
        fasta = self._mock_fasta(tmp_path)
        mock_snakemake.return_value = MagicMock(success=True)

        self.ra.run(sequences_fasta=fasta, output_file=tmp_path / "results.tsv")

        config = mock_snakemake.call_args[1]["config"]
        assert Path(config["sequences_fasta"]).is_absolute()
