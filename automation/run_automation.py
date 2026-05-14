"""Head-container entrypoint for the ONT basecalling automation.

Invoked as ``python -m automation.run_automation`` by the Fargate Batch head
job that the ``startOntBasecall`` Lambda submits when a
``supplemental/barcodes.tsv`` is uploaded to a delivery prefix.

The wrapper:

1. Parses required CLI args supplied by the Lambda via Batch
   ``containerOverrides.command``.
2. Runs ``nextflow run main.nf -profile batch`` against the GPU queue.
3. On success, calls ``seq_import.samplesheet.generate_samplesheet`` in-process
   to write the samplesheet for downstream ``mgs-workflow`` ingestion.

Nextflow runs with ``check=True``; any failure (subprocess or in-process)
propagates a non-zero exit, which surfaces as a ``FAILED`` job in Batch with
the traceback in CloudWatch.
"""

import argparse
import logging
import subprocess

import boto3
from botocore.config import Config
from seq_import.samplesheet import generate_samplesheet

log = logging.getLogger(__name__)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse the entrypoint's CLI args."""
    parser = argparse.ArgumentParser(
        prog="python -m automation.run_automation",
        description="Run basecall-workflow on a delivery, then emit its samplesheet.",
    )
    parser.add_argument(
        "--delivery", required=True,
        help="Delivery name, e.g. NAO-ONT-YYYYMMDD-LIBRARY",
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

    Supplies every param that ``configs/basecall.config`` leaves for the caller
    to provide (the commented-out lines tagged ``// fill ... and uncomment``).
    ``duplex`` and ``demux`` are left to their config defaults (``false`` /
    ``true``). ``-profile batch`` engages the AWS Batch executor + Fusion FS
    settings in ``configs/profiles.config``.
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
    ]


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

    log.info("Generating samplesheet for delivery %s in bucket %s", args.delivery, args.base_bucket)
    s3_client = boto3.client("s3", config=Config(max_pool_connections=50))
    output_path = generate_samplesheet(s3_client, args.delivery, bucket=args.base_bucket)
    log.info("Samplesheet written to: %s", output_path)


if __name__ == "__main__":
    main()
