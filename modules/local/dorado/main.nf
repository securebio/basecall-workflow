// Basecall Nanopore pod5 files
process BASECALL_POD_5_SIMPLEX {
    label "dorado"
    label "basecall"
    accelerator 1

    input:
        tuple path(pod5), val(division)
        val kit
        val nanopore_run

    output:
        tuple path("*.bam"), val(division), emit: bam
        tuple path("sequencing_summary_*.txt"), val(division), emit: summary

    shell:
        '''
        nanopore_run=!{nanopore_run}

        # Dorado basecalling
        dorado basecaller sup !{pod5} --kit-name !{kit} > ${nanopore_run}-!{division}.bam

        dorado summary ${nanopore_run}-!{division}.bam > sequencing_summary_${nanopore_run}-!{division}.txt
        '''
}

process BASECALL_POD_5_DUPLEX {
    label "dorado"
    label "basecall"
    accelerator 1

    input:
        tuple path(pod5), val(division)
        val kit
        val nanopore_run

    output:
        tuple path("*.bam"), val(division), emit: bam
        tuple path("sequencing_summary_*.txt"), val(division), emit: summary

    shell:
        '''
        nanopore_run=!{nanopore_run}

        # Dorado basecalling
        dorado duplex sup !{pod5} > ${nanopore_run}-!{division}.bam

        dorado summary ${nanopore_run}-!{division}.bam > sequencing_summary_${nanopore_run}-!{division}.txt
        '''
}

// Demultiplex basecalled BAM files
process DEMUX_POD_5 {
    label "dorado"
    label "demux"
    accelerator 1

    input:
        tuple path(bam), val(division)
        val kit
        val nanopore_run
        val valid_barcodes
    output:
        // Parent directory is the merge key for downstream groupTuple:
        // per-barcode dirs for valid barcodes, plus an `unclassified/` dir
        // that absorbs both true-unclassified and faulty-barcode reads.
        path('demux_out/*/*.bam'), emit: demux_bam, optional: true

    shell:
        '''
        nanopore_run=!{nanopore_run}
        division=!{division}
        # Store barcodes in a properly quoted variable
        barcodes="!{valid_barcodes}"

        # Turn the barcodes into a proper array by removing brackets and splitting on comma
        barcodes_array=($(echo "$barcodes" | tr -d '[]' | tr ',' ' '))

        mkdir -p tmp_demux

        # Demultiplex
        dorado demux --no-classify --output-dir tmp_demux/ !{bam}

        if [ "$(ls -A tmp_demux/)" ]; then
            for f in tmp_demux/*.bam; do
                demux_id=$(basename "$f" .bam | awk -F '_' '{print $NF}')
                demux_id=${demux_id#barcode}

                if [[ " ${barcodes_array[@]} " =~ " ${demux_id} " ]]; then
                    echo "Processing file: $f with Demux ID: ${demux_id}"
                    merge_key="${demux_id}"
                    new_name="${nanopore_run}-${demux_id}-${division}.bam"
                elif [[ "$f" == *"unclassified"* ]]; then
                    echo "Processing unclassified file: $f"
                    merge_key="unclassified"
                    new_name="${nanopore_run}-unclassified-${division}.bam"
                else
                    echo "Processing wrong barcode: $f"
                    merge_key="unclassified"
                    new_name="${nanopore_run}-faulty-barcode-${demux_id}-${division}.bam"
                fi

                mkdir -p "demux_out/${merge_key}"
                mv "$f" "demux_out/${merge_key}/${new_name}"
            done
        else
            echo "No files to process in tmp_demux/"
        fi

        rmdir tmp_demux 2>/dev/null || true
        '''
}
