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
    Int disk_gb = ceil(1.5*(size(bam, "GiB")))
    String mem = "16 GB"
    Int threads = 16
    Int cpu = (threads)

    command <<<
        echo "Input BAM: ~{bam}"
        echo "Input BED: ~{bed}"

        mkdir -p split_bams split_beds

        # Split BED
        sed 's/^MT/M/' ~{bed} > tmp.bed && mv tmp.bed ~{bed}
        echo "splitting BED for ~{fam_member}"
        awk -v outdir="split_beds" '{ print > (outdir "/" $1 ".bed") }' "~{bed}"

        echo "BEDs created:"
        ls split_beds

        # Split BAM
        for chr in {1..22} M X; do
            echo "Extracting BAM for $chr"
            samtools view -@ 16 -b "~{bam}" "$chr" > "split_bams/$chr.bam"
            samtools index -@ 16 "split_bams/$chr.bam" 
        done

        echo "BAMs created: "
        ls split_bams

    >>>

    output {
        Array[File] bam_array = glob("split_bams/*.bam")
        Array[File] bai_array = glob("split_bams/*.bai")
        Array[File] bed_array = glob("split_beds/*.bed")
    }

    runtime {
        cpu: cpu
        gpu: false
        memory: "${mem}"
        disks: "local-disk ${disk_gb} SSD"
        docker: "${dockerSamtools}"
    }
}
