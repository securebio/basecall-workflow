"""Head-container entrypoint for the ONT basecalling automation.

Invoked as `python -m automation.run_automation` by the Fargate Batch head job
that the `startOntBasecall` Lambda submits when `supplemental/barcodes.tsv` is
uploaded to a delivery prefix.
"""

import logging
import os
import re
import subprocess
import sys

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


def load_env():
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


def build_nextflow_cmd(delivery, kit, aws_queue, base_bucket, work_bucket):
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


def build_samplesheet_cmd(delivery):
    return ["python", "-m", "seq_import", "samplesheet", "--delivery", delivery]


def main():
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
