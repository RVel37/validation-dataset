version 1.0


task split_samples {
    input {
        File bam
        File bai
        File bed
        String fam_member
    }

    # dynamic instance
    Int disk_gb = ceil( 2* (size(bam, "GiB"))) + 2
    String mem = "8 GB"
    Int threads = 4
    Int cpu = (threads)/2

    command <<<

        sample=~{fam_member}

        mkdir -p split_beds/$sample
        mkdir -p split_bams/$sample

        # Split BED by chromosome
        awk -v outdir="split_beds/$sample" '{ print > (outdir "/" $1 ".bed") }' "~{bed}"

        # Split BAM by chromosome + index
        for chrom in $(samtools idxstats "~{bam}" | cut -f1 | grep -v '\*' | sort -V); do
            bam_out="split_bams/$sample/${chrom}.bam"
            bed_out="split_beds/$sample/${chrom}.bed"
            samtools view -b "~{bam}" "$chrom" > "$bam_out"
            samtools index "$bam_out"
        done

    >>>

    output {
        Array[File] bam_array = sort(glob("split_bams/~{fam_member}/*.bam"))
        Array[File] bed_array = sort(glob("split_beds/~{fam_member}/*.bed"))
    }

    runtime {
        cpu: cpu
        gpu: false
        memory: "${mem}"
        disks: "local-disk ${disk_gb} SSD"
    }
}
