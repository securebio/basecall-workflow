"""Head-container entrypoint for the ONT basecalling automation.

Invoked as ``python -m automation.run_automation`` by the Fargate Batch head
job that the ``startOntBasecall`` Lambda submits when a
``supplemental/barcodes.tsv`` is uploaded to a delivery prefix.

The wrapper:

1. Parses required CLI args supplied by the Lambda via Batch
   ``containerOverrides.command``.
2. Validates ``--delivery`` against a conservative regex (defense in depth —
   the Lambda also validates).
3. Runs ``nextflow run main.nf -profile batch`` against the GPU queue.
4. On success, runs ``python -m seq_import samplesheet --delivery <delivery>``
   to write the samplesheet for downstream ``mgs-workflow`` ingestion.

Both subprocesses use ``check=True``; any failure propagates a non-zero exit,
which surfaces as a ``FAILED`` job in Batch with the traceback in CloudWatch.
"""

import argparse
import logging
import re
import subprocess

DUPLEX = "false"
DEMUX = "true"

DELIVERY_RE = re.compile(r"^[A-Za-z0-9_-]+$")

log = logging.getLogger(__name__)


def delivery_type(value: str) -> str:
    """argparse type validator for ``--delivery``."""
    if not DELIVERY_RE.match(value):
        raise argparse.ArgumentTypeError(
            f"{value!r} does not match {DELIVERY_RE.pattern}"
        )
    return value


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse the entrypoint's CLI args."""
    parser = argparse.ArgumentParser(
        prog="python -m automation.run_automation",
        description="Run basecall-workflow on a delivery, then emit its samplesheet.",
    )
    parser.add_argument(
        "--delivery", type=delivery_type, required=True,
        help="Delivery name, e.g. NAO-ONT-YYYYMMDD-LIBRARY (must match [A-Za-z0-9_-]+)",
    )
    parser.add_argument(
        "--kit", required=True,
        help="ONT kit name, e.g. SQK-RPB114-24",
    )
    parser.add_argument(
        "--aws-queue", required=True,
        help="AWS Batch GPU queue for child basecalling jobs",
    )
    parser.add_argument(
        "--base-bucket", required=True,
        help="S3 bucket holding the delivery (raw/, supplemental/, metadata/)",
    )
    parser.add_argument(
        "--work-bucket", required=True,
        help="S3 bucket for Nextflow's working directory",
    )
    return parser.parse_args(argv)


def build_nextflow_cmd(
    delivery: str,
    kit: str,
    aws_queue: str,
    base_bucket: str,
    work_bucket: str,
) -> list[str]:
    """Build the argv for ``nextflow run main.nf``.

    Supplies every param left commented-out in ``configs/basecall.config`` as a
    ``--<param>`` flag. ``-profile batch`` engages the AWS Batch executor +
    Fusion FS settings in ``configs/profiles.config``.
    """
    return [
        "nextflow", "run", "/workflow/main.nf",
        "-profile", "batch",
        "--nanopore_run", delivery,
        "--kit", kit,
        "--aws_queue", aws_queue,
        "--base_dir", f"s3://{base_bucket}/{delivery}",
        "--work_dir", f"s3://{work_bucket}/{delivery}",
        "--barcodes", f"s3://{base_bucket}/{delivery}/supplemental/barcodes.tsv",
        "--duplex", DUPLEX,
        "--demux", DEMUX,
    ]


def build_samplesheet_cmd(delivery: str) -> list[str]:
    """Build the argv for the ``seq_import samplesheet`` CLI.

    ``seq_import`` defaults to ``--bucket nao-restricted``, which matches our
    ``--base-bucket``, so we don't pass ``--bucket`` explicitly — let
    seq_import own that contract.
    """
    return ["python", "-m", "seq_import", "samplesheet", "--delivery", delivery]


def main() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )
    args = parse_args()

    nextflow_cmd = build_nextflow_cmd(
        args.delivery,
        args.kit,
        args.aws_queue,
        args.base_bucket,
        args.work_bucket,
    )
    log.info("Running nextflow: %s", " ".join(nextflow_cmd))
    subprocess.run(nextflow_cmd, check=True, cwd="/workflow")

    samplesheet_cmd = build_samplesheet_cmd(args.delivery)
    log.info("Generating samplesheet: %s", " ".join(samplesheet_cmd))
    subprocess.run(samplesheet_cmd, check=True)


if __name__ == "__main__":
    main()
