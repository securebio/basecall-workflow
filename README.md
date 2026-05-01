# Nanopore Basecalling Workflow

This Nextflow pipeline is designed to process Oxford Nanopore raw signal data (POD5 files) through basecalling and optional demultiplexing steps. It supports both simplex and duplex basecalling modes using Dorado.

## Pipeline Description

### Overview

The pipeline consists of a single workflow that processes Nanopore POD5 files through several phases:

1. A **basecalling phase** using Dorado in either simplex or duplex mode
2. An optional **demultiplexing phase** for barcoded samples
3. A final **conversion phase** to generate FASTQ files from BAM output

### Pipeline Outputs

The workflow produces the following outputs:

1. `raw/`: Directory containing the final FASTQ files
2. `unclassified/`: Directory containing unclassified FASTQ files (only relevant for demultiplexing)

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
