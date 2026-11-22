"""Barcode formatting and demultiplexing read-count QC helpers for basecall-workflow."""

from __future__ import annotations

import re
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

# ONT barcodes 01-96 in zero-padded two-digit format (no NB/RLB/barcode prefix).
VALID_BARCODE_RE = re.compile(r"^(?:0[1-9]|[1-8][0-9]|9[0-6])$")


def parse_and_validate_barcodes(lines: Sequence[str]) -> list[str]:
    """Parse and validate lines from ``barcodes.tsv``.

    Raises ``ValueError`` if the list is empty, contains duplicates, or has any
    entry missing a leading zero (e.g. ``1``) or retaining an ONT prefix
    (e.g. ``NB01``, ``RLB01``, ``barcode01``).
    """
    barcodes: list[str] = []
    seen: set[str] = set()
    for line_no, raw in enumerate(lines, start=1):
        token = raw.strip()
        if not token:
            continue
        if not VALID_BARCODE_RE.fullmatch(token):
            raise ValueError(
                f"Invalid barcode {token!r} on line {line_no}: expected zero-padded "
                "two-digit ONT barcode ID '01'-'96' without 'NB'/'RLB'/'barcode' prefix."
            )
        if token in seen:
            raise ValueError(f"Duplicate barcode entry {token!r} on line {line_no}.")
        seen.add(token)
        barcodes.append(token)

    if not barcodes:
        raise ValueError(
            "Barcode list is empty; expected at least one valid barcode ('01'-'96')."
        )
    return barcodes


def validate_barcodes_file(path: Path) -> list[str]:
    """Read and validate a ``barcodes.tsv`` file from disk."""
    lines = path.read_text(encoding="utf-8").splitlines()
    return parse_and_validate_barcodes(lines)


def evaluate_demux_read_counts(
    expected_barcodes: Sequence[str],
    observed_counts: Mapping[str, int],
    unclassified_reads: int = 0,
    min_included_reads: int = 10,
    max_excluded_fraction: float = 0.05,
) -> dict[str, Any]:
    """Check demultiplexed per-barcode read counts for mis-specified barcodes.tsv.

    Flags two silent-failure conditions:
    1. Any expected barcode with fewer than ``min_included_reads`` reads (including 0).
    2. Any excluded (unexpected) barcode exceeding ``max_excluded_fraction`` of
       total classified reads (or total excluded reads exceeding that threshold).
    """
    if min_included_reads < 0:
        raise ValueError("min_included_reads must be non-negative.")
    if not (0.0 <= max_excluded_fraction <= 1.0):
        raise ValueError("max_excluded_fraction must be between 0.0 and 1.0.")

    expected_clean = parse_and_validate_barcodes(list(expected_barcodes))
    expected_set = set(expected_clean)

    included_counts: dict[str, int] = {bc: 0 for bc in expected_clean}
    excluded_counts: dict[str, int] = {}

    for raw_bc, raw_count in observed_counts.items():
        bc = raw_bc.strip()
        count = int(raw_count)
        if count < 0:
            raise ValueError(f"Negative read count ({count}) for barcode {bc!r}.")
        if bc in expected_set:
            included_counts[bc] = included_counts.get(bc, 0) + count
        elif count > 0:
            excluded_counts[bc] = excluded_counts.get(bc, 0) + count

    total_included = sum(included_counts.values())
    total_excluded = sum(excluded_counts.values())
    total_classified = total_included + total_excluded

    failures: list[str] = []
    if total_classified == 0:
        failures.append(
            f"Zero classified reads observed across all expected barcodes {expected_clean} "
            f"(unclassified={unclassified_reads})."
        )
    else:
        for bc in expected_clean:
            cnt = included_counts[bc]
            if cnt < min_included_reads:
                failures.append(
                    f"Included barcode {bc!r} has only {cnt} reads (< min_included_reads={min_included_reads})."
                )
        for bc, cnt in sorted(excluded_counts.items()):
            frac = cnt / total_classified
            if frac > max_excluded_fraction:
                failures.append(
                    f"Excluded barcode {bc!r} accounts for {cnt}/{total_classified} "
                    f"({frac:.1%}) of classified reads (> max_excluded_fraction={max_excluded_fraction:.1%})."
                )
        excluded_frac = total_excluded / total_classified
        if (
            excluded_counts
            and excluded_frac > max_excluded_fraction
            and not any("Excluded barcode" in msg for msg in failures)
        ):
            failures.append(
                f"Combined excluded barcodes account for {total_excluded}/{total_classified} "
                f"({excluded_frac:.1%}) of classified reads (> max_excluded_fraction={max_excluded_fraction:.1%})."
            )

    return {
        "passed": len(failures) == 0,
        "failures": failures,
        "included_counts": included_counts,
        "excluded_counts": excluded_counts,
        "total_included_reads": total_included,
        "total_excluded_reads": total_excluded,
        "total_classified_reads": total_classified,
        "unclassified_reads": int(unclassified_reads),
    }


def check_demux_qc_or_raise(
    expected_barcodes: Sequence[str],
    observed_counts: Mapping[str, int],
    unclassified_reads: int = 0,
    min_included_reads: int = 10,
    max_excluded_fraction: float = 0.05,
) -> dict[str, Any]:
    """Evaluate demux barcode QC and raise ``ValueError`` with details on failure."""
    report = evaluate_demux_read_counts(
        expected_barcodes=expected_barcodes,
        observed_counts=observed_counts,
        unclassified_reads=unclassified_reads,
        min_included_reads=min_included_reads,
        max_excluded_fraction=max_excluded_fraction,
    )
    if not report["passed"]:
        raise ValueError("Demux barcode QC failed: " + "; ".join(report["failures"]))
    return report
