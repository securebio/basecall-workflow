// Basecall Nanopore pod5 files
process BASECALL_POD_5_SIMPLEX {
    label "dorado"
    label "basecall"
    accelerator 1
    memory '16 GB'

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

// Option with no trimming -- can be useful for R&D
process BASECALL_POD_5_NOTRIM {
    label "dorado"
    label "basecall"
    accelerator 1
    memory '16 GB'

    input:
        tuple path(pod5), val(division)
        val nanopore_run
        
    output:
        tuple path("*.bam"), val(division), emit: bam
        tuple path("sequencing_summary_*.txt"), val(division), emit: summary

    shell:
        '''
        nanopore_run=!{nanopore_run}

        # Dorado basecalling
        dorado basecaller sup !{pod5} --no-trim > ${nanopore_run}-!{division}.bam

        dorado summary ${nanopore_run}-!{division}.bam > sequencing_summary_${nanopore_run}-!{division}.txt
        '''
}


process BASECALL_POD_5_DUPLEX {
    label "dorado"
    label "basecall"
    accelerator 1
    memory '16 GB'

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
    memory '16 GB'

    input:
        tuple path(bam), val(division)
        val kit
        val nanopore_run
        val valid_barcodes
    output:
        path('demultiplexed/*'), emit: demux_bam, optional: true
        path('unclassified/*'), emit: unclassified_bam, optional: true

    shell:
        '''
        nanopore_run=!{nanopore_run}
        division=!{division}
        # Store barcodes in a properly quoted variable
        barcodes="!{valid_barcodes}"

        # Turn the barcodes into a proper array by removing brackets and splitting on comma
        barcodes_array=($(echo "$barcodes" | tr -d '[]' | tr ',' ' '))
        
        # Create unclassified and demultiplexed dirs
        mkdir -p unclassified
        mkdir -p demultiplexed

        # Demultiplex
        dorado demux --no-classify --output-dir demultiplexed/ !{bam}

        # Rename output files
        if [ "$(ls -A demultiplexed/)" ]; then
            for f in demultiplexed/*; do
                # Extract demux_id from filename
                demux_id=$(basename "$f" .bam | awk -F '_' '{print $NF}')
                demux_id=${demux_id#barcode}

                # Check if demux_id is in valid_barcodes
                if [[ " ${barcodes_array[@]} " =~ " ${demux_id} " ]]; then
                    echo "Processing file: $f with Demux ID: ${demux_id}"
                    mv "$f" "demultiplexed/${nanopore_run}-${demux_id}-${division}.bam"
                elif [[ "$f" == *"unclassified"* ]]; then
                    echo "Processing unclassified file: $f"
                    mv "$f" "unclassified/${nanopore_run}-unclassified-${division}.bam"
                else
                    echo "Processing wrong barcode: $f"
                    mv "$f" "unclassified/${nanopore_run}-faulty-barcode-${demux_id}-${division}.bam"
                fi
            done
        else
            echo "No files to process in demultiplexed/"
        fi
        '''
}
