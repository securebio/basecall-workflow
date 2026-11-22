"""Unit and CLI tests for barcode formatting validation and demux read-count QC."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from automation.barcode_qc import (
    check_demux_qc_or_raise,
    evaluate_demux_read_counts,
    parse_and_validate_barcodes,
    validate_barcodes_file,
)
from bin import validate_barcodes as validate_cli


class TestBarcodeFormatValidation:
    def test_accepts_zero_padded_ont_barcodes(self) -> None:
        assert parse_and_validate_barcodes(["01", " 02 ", "", "12", "96"]) == [
            "01",
            "02",
            "12",
            "96",
        ]

    @pytest.mark.parametrize(
        "bad_entry",
        ["1", "9", "00", "97", "NB01", "RLB04", "barcode07", "bc01", "abc"],
    )
    def test_rejects_unpadded_or_prefixed_barcodes(self, bad_entry: str) -> None:
        with pytest.raises(ValueError, match="zero-padded"):
            parse_and_validate_barcodes(["01", bad_entry])

    def test_rejects_empty_and_duplicate_entries(self) -> None:
        with pytest.raises(ValueError, match="empty"):
            parse_and_validate_barcodes(["", "   "])
        with pytest.raises(ValueError, match="Duplicate"):
            parse_and_validate_barcodes(["01", "02", "01"])

    def test_validate_barcodes_file(self, tmp_path: Path) -> None:
        tsv = tmp_path / "barcodes.tsv"
        tsv.write_text("01\n02\n\n03\n", encoding="utf-8")
        assert validate_barcodes_file(tsv) == ["01", "02", "03"]


class TestDemuxReadCountQC:
    def test_clean_demux_passes(self) -> None:
        rep = check_demux_qc_or_raise(
            expected_barcodes=["01", "02"],
            observed_counts={"01": 500, "02": 480, "03": 5},
            unclassified_reads=120,
            min_included_reads=10,
            max_excluded_fraction=0.05,
        )
        assert rep["passed"] is True
        assert rep["total_classified_reads"] == 985

    def test_fails_when_included_barcode_missing_or_too_small(self) -> None:
        rep = evaluate_demux_read_counts(
            expected_barcodes=["01", "02"],
            observed_counts={"01": 500, "02": 3},
            min_included_reads=10,
        )
        assert rep["passed"] is False
        assert any("02" in msg and "3 reads" in msg for msg in rep["failures"])

        with pytest.raises(ValueError, match="02"):
            check_demux_qc_or_raise(
                expected_barcodes=["01", "02"],
                observed_counts={"01": 500},
                min_included_reads=10,
            )

    def test_fails_when_excluded_barcode_has_large_share(self) -> None:
        rep = evaluate_demux_read_counts(
            expected_barcodes=["01", "02"],
            observed_counts={"01": 400, "02": 400, "07": 200},
            min_included_reads=10,
            max_excluded_fraction=0.05,
        )
        assert rep["passed"] is False
        assert any("Excluded barcode '07'" in msg for msg in rep["failures"])

    def test_fails_when_all_reads_unclassified_or_zero_classified(self) -> None:
        with pytest.raises(ValueError, match="Zero classified reads"):
            check_demux_qc_or_raise(
                expected_barcodes=["01", "02"],
                observed_counts={},
                unclassified_reads=1000,
            )


class TestBarcodeQcCli:
    def test_cli_validate_and_check_counts(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        good_tsv = tmp_path / "good.tsv"
        good_tsv.write_text("01\n02\n", encoding="utf-8")
        assert validate_cli.main(["validate-barcodes", str(good_tsv)]) == 0
        assert json.loads(capsys.readouterr().out)["barcodes"] == ["01", "02"]

        bad_tsv = tmp_path / "bad.tsv"
        bad_tsv.write_text("1\nNB02\n", encoding="utf-8")
        assert validate_cli.main(["validate-barcodes", str(bad_tsv)]) == 1
        assert "ERROR:" in capsys.readouterr().err

        counts_tsv = tmp_path / "counts.tsv"
        counts_tsv.write_text(
            "barcode\treads\n01\t200\n02\t200\n05\t150\n", encoding="utf-8"
        )
        assert (
            validate_cli.main(
                [
                    "check-demux-counts",
                    "--barcodes",
                    str(good_tsv),
                    "--counts",
                    str(counts_tsv),
                ]
            )
            == 1
        )
        assert "Excluded barcode '05'" in capsys.readouterr().err
