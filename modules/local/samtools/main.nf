process BAM_TO_FASTQ {
    label "samtools"

    input:
        path(bam)
    output:
        path '*.fastq.gz'

    shell:
        '''
        # Get base name without .bam extension
        base_name=$(basename !{bam} .bam)

        # Convert BAM to FASTQ and gzip
        samtools fastq !{bam} | gzip -c > "${base_name}.fastq.gz"
        '''
}

process MERGE_BAMS {
    label "samtools"

    input:
        tuple val(out_name), path(bam_files)
    output:
        path("${out_name}.bam")
    shell:
        '''
        samtools merge -r -o !{out_name}.bam !{bam_files}
        '''
}

process CHECK_DEMUX_BARCODES {
    label "samtools"

    input:
        path(demux_bams)
        val(expected_barcodes)
        val(nanopore_run)
        val(min_included_reads)
        val(max_excluded_fraction)
        val(fail_on_qc)
    output:
        path("demux_barcode_qc.tsv"), emit: qc_report

    shell:
        '''
        set -euo pipefail
        nanopore_run="!{nanopore_run}"
        expected_csv="!{expected_barcodes.join(',')}"
        min_included="!{min_included_reads}"
        max_excluded_frac="!{max_excluded_fraction}"
        fail_on_qc="!{fail_on_qc}"

        raw_counts="raw_barcode_counts.tsv"
        : > "$raw_counts"

        for bam in !{demux_bams}; do
            [ -e "$bam" ] || continue
            bname=$(basename "$bam" .bam)
            reads=$(samtools view -c "$bam")
            prefix="${nanopore_run}-"
            rest="${bname#"$prefix"}"
            if [[ "$rest" == faulty-barcode-* ]]; then
                sub="${rest#faulty-barcode-}"
                bc="${sub%-*}"
                printf "excluded\t%s\t%s\n" "${bc}" "${reads}" >> "$raw_counts"
            elif [[ "$rest" == unclassified-* ]]; then
                printf "unclassified\tunclassified\t%s\n" "${reads}" >> "$raw_counts"
            else
                bc="${rest%-*}"
                printf "included\t%s\t%s\n" "${bc}" "${reads}" >> "$raw_counts"
            fi
        done

        awk -F'\t' -v OFS='\t' \
            -v expected_csv="$expected_csv" \
            -v min_inc="$min_included" \
            -v max_exc_frac="$max_excluded_frac" \
            '
            BEGIN {
                n_exp = split(expected_csv, exp_arr, ",");
                for (i = 1; i <= n_exp; i++) {
                    if (exp_arr[i] != "") {
                        exp_seen[exp_arr[i]] = 1;
                        inc_cnt[exp_arr[i]] = 0;
                    }
                }
                total_inc = 0;
                total_exc = 0;
                unclass = 0;
            }
            {
                kind = $1; bc = $2; n = $3 + 0;
                if (kind == "included") {
                    inc_cnt[bc] += n;
                    total_inc += n;
                } else if (kind == "excluded") {
                    exc_cnt[bc] += n;
                    total_exc += n;
                } else if (kind == "unclassified") {
                    unclass += n;
                }
            }
            END {
                total_cls = total_inc + total_exc;
                print "category", "barcode", "reads", "fraction_of_classified", "status";
                failed = 0;
                if (total_cls == 0) {
                    print "summary", "ALL", 0, "0.0000", "FAIL_ZERO_CLASSIFIED";
                    failed = 1;
                }
                for (i = 1; i <= n_exp; i++) {
                    bc = exp_arr[i];
                    if (bc == "") continue;
                    c = inc_cnt[bc] + 0;
                    frac = (total_cls > 0 ? c / total_cls : 0.0);
                    st = (c < min_inc ? "FAIL_LOW_INCLUDED" : "OK");
                    if (st != "OK") failed = 1;
                    printf "included\t%s\t%d\t%.4f\t%s\n", bc, c, frac, st;
                }
                for (bc in exc_cnt) {
                    c = exc_cnt[bc] + 0;
                    frac = (total_cls > 0 ? c / total_cls : 1.0);
                    st = (frac > max_exc_frac ? "FAIL_HIGH_EXCLUDED" : "OK");
                    if (st != "OK") failed = 1;
                    printf "excluded\t%s\t%d\t%.4f\t%s\n", bc, c, frac, st;
                }
                printf "unclassified\tunclassified\t%d\tNA\tINFO\n", unclass;
                if (failed) {
                    exit 42;
                }
            }
            ' "$raw_counts" > demux_barcode_qc.tsv || rc=$?

        rc=${rc:-0}
        if [ "$rc" -ne 0 ]; then
            cat demux_barcode_qc.tsv >&2
            if [ "$fail_on_qc" = "true" ]; then
                echo "ERROR: Demultiplexing barcode QC failed (possible incorrect barcodes.tsv)." >&2
                exit 1
            fi
        fi
        '''
}
