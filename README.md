# Nanopore Basecalling Workflow

This Nextflow pipeline is designed to process Oxford Nanopore raw signal data (POD5 files) through basecalling and optional demultiplexing steps. It supports both simplex and duplex basecalling modes using Dorado.

## Pipeline Description

### Overview

The pipeline consists of a single workflow that processes Nanopore POD5 files through several phases:

1. A **basecalling phase** using Dorado in either simplex or duplex mode
2. An optional **demultiplexing phase** for barcoded samples
3. A final **conversion phase** to generate FASTQ files from BAM output

### Pipeline Outputs

All output FASTQs are published to a single `raw/` directory.

When `demux = true`, the directory contains one FASTQ per barcode listed in `barcodes.txt`, plus one for unclassified reads:

- `${nanopore_run}-${barcode}_SE.fastq.gz` (e.g., `${nanopore_run}-01_SE.fastq.gz`)
- `${nanopore_run}-unclassified_SE.fastq.gz` (also includes reads that dorado tagged with a barcode not in `barcodes.txt`)

The `_SE` suffix denotes single-end reads, distinguishing these files from the `_R1` / `_R2` suffix convention used for paired-end reads.

## Using the Workflow

### Installation & Setup

1. Install Nextflow (23.04.0+)
2. Install Docker
3. Set up [AWS BATCH](https://github.com/naobservatory/mgs-workflow/tree/master#:~:text=The%20batch%20profile%20is,your%20Batch%20job%20queue.)
4. Clone this repository

### Running the Pipeline

Basic usage:

Create a new directory, name it after the delivery, copy in basecall.config as nextflow.config, and set the parameters. Params:

- duplex
  - Duplex basecalling or no? You can't combine duplex and demux
- demux
  - Demultiplex basecalling output?
- nanopore_run
  - Name of run/delivery
- kit
  - Name of ONT kit, needed for demux'ing

Additionally, add a barcodes.txt file to the directory, containing the barcodes to be demultiplexed, in the format:

```
01
02
12
...
```

Once that is done, you can switch into the directory and run

```bash
nextflow run .. -resume
```

## Automation

The `automation/` directory contains the head container that wraps this workflow for automated execution as part of the ONT basecalling automation pipeline in [`det-terraform-production`](https://github.com/securebio/det-terraform-production). The container bundles Nextflow, the workflow source, and `seq_import` (the upload helper from [`nao-mgs-import`](https://github.com/securebio/nao-mgs-import)).

- `automation/Dockerfile` — head container image
- `automation/environment.yml` — micromamba env spec (Python, Nextflow, AWS CLI, git)
- `.github/workflows/ecr-push.yml` — publishes the image to ECR on push to `main`
- `.github/workflows/docker-build.yml` — PR-time build smoke test

### Container entrypoint

The image runs `python -m automation.run_automation`, which invokes `nextflow run main.nf` against the AWS Batch GPU queue and, on success, runs `python -m seq_import samplesheet --delivery $DELIVERY` to write the samplesheet to `s3://$BASE_BUCKET/$DELIVERY/metadata/samplesheet.csv`.

Required environment variables (passed by the `startOntBasecall` Lambda via Batch `containerOverrides`):

- `DELIVERY` — delivery name, e.g. `NAO-ONT-YYYYMMDD-LIBRARY` (must match `[A-Za-z0-9_-]+`)
- `KIT` — ONT kit name, e.g. `SQK-RPB114-24`
- `AWS_QUEUE` — AWS Batch GPU queue for child basecalling jobs
- `BASE_BUCKET` — S3 bucket holding the delivery (`raw/`, `supplemental/`, `metadata/`)
- `WORK_BUCKET` — S3 bucket for Nextflow's working directory

Local smoke test:

```bash
docker run --rm \
  -e DELIVERY=NAO-ONT-YYYYMMDD-LIBRARY \
  -e KIT=SQK-RPB114-24 \
  -e AWS_QUEUE=<gpu-queue> \
  -e BASE_BUCKET=nao-restricted \
  -e WORK_BUCKET=sb-det-ont-basecall-work \
  basecall-workflow
```

### Build-time authentication

Because the image pip-installs `seq_import` from the private `nao-mgs-import` repo, the build needs a GitHub token. CI mints a short-lived one via the `sbd-mgs-import-reader` GitHub App (App ID in `vars.IMPORT_READER_APP_ID`, private key in `secrets.IMPORT_READER_PRIVATE_KEY`) and passes it to `docker build` as a BuildKit secret, so it never lands in image layers. To build locally, supply any token with read access to `nao-mgs-import`:

```bash
GH_TOKEN=$(gh auth token) docker build \
  --secret id=gh_token,env=GH_TOKEN \
  -f automation/Dockerfile -t basecall-workflow .
```
