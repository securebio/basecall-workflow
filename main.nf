/***********************************************************************************************
| WORKFLOW: BASECALLING NANOPORE SQUIGGLE DATA |
***********************************************************************************************/

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { BASECALL_POD_5_SIMPLEX } from "./modules/local/dorado"
include { BASECALL_POD_5_DUPLEX } from "./modules/local/dorado"
include { DEMUX_POD_5 } from "./modules/local/dorado"
include { BAM_TO_FASTQ } from "./modules/local/samtools"
include { MERGE_BAMS } from "./modules/local/samtools"
include { CHECK_DEMUX_BARCODES } from "./modules/local/samtools"

/*****************
| MAIN WORKFLOWS |
*****************/

// Complete primary workflow
workflow {
    main:
    // Start time
    start_time = new Date()
    start_time_str = start_time.format("YYYY-MM-dd HH:mm:ss z (Z)")

    // Batching
    pod5_ch = channel.fromPath("${params.base_dir}/pod5/*")

    // file -> tuple(file, division)
    pod5_ch = pod5_ch.collect(flat: false, sort: true)
        .flatMap { files ->
        files.withIndex().collect { file, index ->
            tuple(file, String.format("div%04d", index + 1))
        }
    }

    // Barcodes (fail fast before GPU basecalling if formatting is invalid)
    def raw_barcodes = file(params.barcodes).readLines()
        .collect { it.trim() }
        .findAll { it }
    if (params.demux && !params.duplex) {
        if (raw_barcodes.isEmpty()) {
            error("Barcode list '${params.barcodes}' is empty.")
        }
        def seen_barcodes = [] as Set
        raw_barcodes.each { bc ->
            if (!(bc ==~ /^(?:0[1-9]|[1-8][0-9]|9[0-6])$/)) {
                error("Invalid barcode '${bc}' in '${params.barcodes}': expected zero-padded two-digit ONT ID (e.g. '01', without 'NB'/'RLB'/'barcode' prefix).")
            }
            if (!seen_barcodes.add(bc)) {
                error("Duplicate barcode '${bc}' in '${params.barcodes}'.")
            }
        }
    }
    barcodes_ch = raw_barcodes

    // Basecalling
    if (params.duplex) {
        bam_ch = BASECALL_POD_5_DUPLEX(pod5_ch, params.nanopore_run)
        final_bam_ch = bam_ch.bam.flatten()
    } else {
        bam_ch = BASECALL_POD_5_SIMPLEX(pod5_ch, params.kit, params.nanopore_run)
        if (params.demux) {
            demux_ch = DEMUX_POD_5(bam_ch.bam, params.nanopore_run, barcodes_ch)

            CHECK_DEMUX_BARCODES(
                demux_ch.demux_bam.flatten().collect(),
                barcodes_ch,
                params.nanopore_run,
                params.min_included_barcode_reads,
                params.max_excluded_barcode_fraction,
                params.fail_on_barcode_qc
            )

            // Each BAM's parent directory name is its merge key, set by DEMUX_POD_5.
            merge_input_ch = demux_ch.demux_bam.flatten()
                .map { bam -> tuple("${params.nanopore_run}-${bam.parent.name}_SE", bam) }
                .groupTuple()

            final_bam_ch = MERGE_BAMS(merge_input_ch)
        }
    }

    // Convert to FASTQ
    fastq_ch = BAM_TO_FASTQ(final_bam_ch)

    publish:
        fastq_ch = fastq_ch
}


output {
     fastq_ch {
        path "raw"
        tags nextflow_file_class: "publish", "nextflow.io/temporary": "false"
    }
}
