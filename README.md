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
