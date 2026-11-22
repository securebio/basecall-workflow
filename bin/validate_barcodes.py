#!/usr/bin/env python3
"""CLI entrypoint for validating ONT barcodes.tsv and demultiplexing read-count QC."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from automation.barcode_qc import check_demux_qc_or_raise, validate_barcodes_file  # noqa: E402


def _load_counts_tsv(counts_path: Path) -> tuple[dict[str, int], int]:
    observed: dict[str, int] = {}
    unclassified = 0
    for line in counts_path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#") or line.lower().startswith("barcode\t"):
            continue
        parts = line.split("\t")
        if len(parts) < 2:
            raise ValueError(f"Malformed counts TSV row: {line!r}")
        bc, count_str = parts[0].strip(), parts[1].strip()
        count = int(count_str)
        if bc == "unclassified":
            unclassified += count
        else:
            observed[bc] = observed.get(bc, 0) + count
    return observed, unclassified


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Validate basecall-workflow barcodes and demux counts."
    )
    sub = parser.add_subparsers(dest="command", required=True)

    pv = sub.add_parser("validate-barcodes", help="Validate barcodes.tsv formatting.")
    pv.add_argument("barcodes_tsv", type=Path)

    pc = sub.add_parser(
        "check-demux-counts", help="Verify included vs excluded barcode read counts."
    )
    pc.add_argument("--barcodes", type=Path, required=True)
    pc.add_argument("--counts", type=Path, required=True)
    pc.add_argument("--min-included-reads", type=int, default=10)
    pc.add_argument("--max-excluded-fraction", type=float, default=0.05)

    args = parser.parse_args(argv)
    try:
        if args.command == "validate-barcodes":
            barcodes = validate_barcodes_file(args.barcodes_tsv)
            print(json.dumps({"valid": True, "barcodes": barcodes}))
            return 0
        if args.command == "check-demux-counts":
            expected = validate_barcodes_file(args.barcodes)
            observed, unclassified = _load_counts_tsv(args.counts)
            report = check_demux_qc_or_raise(
                expected_barcodes=expected,
                observed_counts=observed,
                unclassified_reads=unclassified,
                min_included_reads=args.min_included_reads,
                max_excluded_fraction=args.max_excluded_fraction,
            )
            print(json.dumps(report))
            return 0
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
