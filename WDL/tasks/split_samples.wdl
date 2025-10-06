version 1.0


task split_samples {
    input {
        File bam
        File bai
        File bed
        String fam_member
        String dockerSamtools
    }

    # dynamic instance
    Int disk_gb = ceil( 2* (size(bam, "GiB"))) + 2
    String mem = "16 GB"
    Int threads = 16
    Int cpu = (threads)/2

    command <<<
        echo "Input BAM: ~{bam}"
        echo "Input BED: ~{bed}"

        mkdir -p split_bams split_beds

        # Split BED
        echo "splitting BED for ~{fam_member}"
        awk -v outdir="split_beds" '{ print > (outdir "/" $1 ".bed") }' "~{bed}"

        # Split BAM
        echo "splitting BAM for ~{fam_member}"
        samtools split -f "split_bams/%!.bam" "~{bam}"
        
        echo "BAMs created: "
        ls split_bams

    >>>

    output {
        Array[File] bam_array = glob("split_bams/~{fam_member}/*.bam")
        Array[File] bed_array = glob("split_beds/~{fam_member}/*.bed")
    }

    runtime {
        cpu: cpu
        gpu: false
        memory: "${mem}"
        disks: "local-disk ${disk_gb} SSD"
        docker: "${dockerSamtools}"
    }
}
