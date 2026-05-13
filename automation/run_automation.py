"""Head-container entrypoint for the ONT basecalling automation.

Invoked as ``python -m automation.run_automation`` by the Fargate Batch head
job that the ``startOntBasecall`` Lambda submits when a
``supplemental/barcodes.tsv`` is uploaded to a delivery prefix.

The wrapper:

1. Reads required env vars set by the Lambda via Batch ``containerOverrides``.
2. Validates ``DELIVERY`` against a conservative regex (defense in depth — the
   Lambda also validates).
3. Runs ``nextflow run main.nf -profile batch`` against the GPU queue.
4. On success, runs ``python -m seq_import samplesheet --delivery $DELIVERY``
   to write the samplesheet for downstream ``mgs-workflow`` ingestion.

Both subprocesses use ``check=True``; any failure propagates a non-zero exit,
which surfaces as a ``FAILED`` job in Batch with the traceback in CloudWatch.
"""

import logging
import os
import re
import subprocess

DUPLEX = "false"
DEMUX = "true"

DELIVERY_RE = re.compile(r"^[A-Za-z0-9_-]+$")

REQUIRED_ENV_VARS = (
    "DELIVERY",
    "KIT",
    "AWS_QUEUE",
    "BASE_BUCKET",
    "WORK_BUCKET",
)

log = logging.getLogger(__name__)


def load_env() -> dict[str, str]:
    """Read and validate required env vars; exit non-zero if any are missing or invalid."""
    missing = [v for v in REQUIRED_ENV_VARS if not os.environ.get(v)]
    if missing:
        raise SystemExit(
            f"Missing required environment variable(s): {', '.join(missing)}"
        )
    delivery = os.environ["DELIVERY"]
    if not DELIVERY_RE.match(delivery):
        raise SystemExit(
            f"DELIVERY {delivery!r} does not match {DELIVERY_RE.pattern}"
        )
    return {v: os.environ[v] for v in REQUIRED_ENV_VARS}


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
    ``BASE_BUCKET``, so we don't pass ``--bucket`` explicitly — let seq_import
    own that contract.
    """
    return ["python", "-m", "seq_import", "samplesheet", "--delivery", delivery]


def main() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )
    env = load_env()

    nextflow_cmd = build_nextflow_cmd(
        env["DELIVERY"],
        env["KIT"],
        env["AWS_QUEUE"],
        env["BASE_BUCKET"],
        env["WORK_BUCKET"],
    )
    log.info("Running nextflow: %s", " ".join(nextflow_cmd))
    subprocess.run(nextflow_cmd, check=True, cwd="/workflow")

    samplesheet_cmd = build_samplesheet_cmd(env["DELIVERY"])
    log.info("Generating samplesheet: %s", " ".join(samplesheet_cmd))
    subprocess.run(samplesheet_cmd, check=True)


if __name__ == "__main__":
    main()
