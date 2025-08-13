/***********************************************************************************************
| WORKFLOW: BASECALLING NANOPORE SQUIGGLE DATA |
***********************************************************************************************/

import groovy.json.JsonOutput
import java.time.LocalDateTime

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { BASECALL_POD_5_SIMPLEX } from "./modules/local/dorado"
include { BASECALL_POD_5_DUPLEX } from "./modules/local/dorado"
include { DEMUX_POD_5 } from "./modules/local/dorado"
include { BAM_TO_FASTQ } from "./modules/local/samtools"
include { MERGE_BAMS } from "./modules/local/samtools"
nextflow.preview.output = true

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

    // Basecalling
    if (params.duplex) {
        bam_ch = BASECALL_POD_5_DUPLEX(pod5_ch, params.kit, params.nanopore_run)
        final_bam_ch = bam_ch.bam.flatten()
    } else {
        if (params.notrim) {
            bam_ch = BASECALL_POD_5_NOTRIM(pod5_ch, params.nanopore_run)
        } else {
            bam_ch = BASECALL_POD_5_SIMPLEX(pod5_ch, params.kit, params.nanopore_run)
        }
        if (params.demux) {
            barcodes_ch = file(params.barcodes).readLines().collect()
            demux_ch = DEMUX_POD_5(bam_ch.bam, params.kit, params.nanopore_run, barcodes_ch)
            classified_bam_ch = demux_ch.demux_bam.flatten()
            unclassified_bam_ch = MERGE_BAMS(demux_ch.unclassified_bam.collect(), params.nanopore_run)
            final_bam_ch = classified_bam_ch.mix(unclassified_bam_ch)
        }
        else {
            final_bam_ch = bam_ch.bam.map { bam, _division -> bam}.flatten()
        }
    }

    // Convert to FASTQ
    fastq_ch = BAM_TO_FASTQ(final_bam_ch, params.nanopore_run)

    publish:
        fastq_ch = fastq_ch
}


output {
     fastq_ch {
        path "raw_trimmed"
        tags nextflow_file_class: "publish", "nextflow.io/temporary": "false"
    }
}
