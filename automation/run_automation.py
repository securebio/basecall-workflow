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
import datetime as dt
import logging
import subprocess
from pathlib import Path

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
    parser.add_argument(
        "--log-bucket", required=True,
        help="S3 bucket for publishing .nextflow.log after the run",
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
        "-c", "/workflow/configs/basecall.config",
        "-profile", "batch",
        "--nanopore_run", delivery,
        "--kit", kit,
        "--aws_queue", aws_queue,
        "--base_dir", f"s3://{base_bucket}/{delivery}",
        "--work_dir", f"s3://{work_bucket}/{delivery}",
        "--barcodes", f"s3://{base_bucket}/{delivery}/supplemental/barcodes.tsv",
    ]


def upload_nextflow_log(
    s3_client, delivery: str, log_bucket: str, log_path: Path = Path("/workflow/.nextflow.log"),
) -> None:
    """Upload .nextflow.log to s3://{log_bucket}/basecall-workflow/automated/{delivery}/{ts}/.

    Logs and swallows all upload errors so it is safe to call from a finally block.
    The error log surfaces the failure in CloudWatch so a misconfigured IAM/bucket gets noticed.
    """
    if not log_path.exists():
        log.warning("No %s to upload", log_path)
        return
    timestamp = dt.datetime.now(dt.UTC).strftime("%Y%m%d_%H%M%S")
    s3_key = f"basecall-workflow/automated/{delivery}/{timestamp}/.nextflow.log"
    try:
        s3_client.upload_file(str(log_path), log_bucket, s3_key)
        log.info("Uploaded nextflow log to s3://%s/%s", log_bucket, s3_key)
    except Exception as e:
        log.exception("Failed to upload %s: %s", log_path, e)


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
    s3_client = boto3.client("s3", config=Config(max_pool_connections=50))
    try:
        subprocess.run(nextflow_cmd, check=True, cwd="/workflow")
    finally:
        # upload_nextflow_log is designed to log-and-swallow all exceptions.
        # (Because an exception in this finally block would prevent 
        # samplesheet generation on a successful Nextflow run, or mask CalledProcessError 
        # from a failed one.)
        upload_nextflow_log(s3_client, args.delivery, args.log_bucket)

    log.info("Generating samplesheet for delivery %s in bucket %s", args.delivery, args.base_bucket)
    output_path = generate_samplesheet(s3_client, args.delivery, bucket=args.base_bucket)
    log.info("Samplesheet written to: %s", output_path)


if __name__ == "__main__":
    main()
