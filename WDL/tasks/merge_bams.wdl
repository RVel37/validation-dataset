version 1.0

task merge_bams {
    input {
        Array[File] bams
        String fam_member
        String dockerSamtools
    }

    Int disk_gb = ceil(1.5 * size(bams, "GiB"))
    String mem = "16 GB"
    Int threads = 8
    Int cpu = 8

    command <<<
        # DEBUG
        echo "List of BAMs:"
        for bam in ~{sep=' ' bams}; do
            echo $bam
        done

        echo "Merging BAMs for ~{fam_member}"
        samtools merge -@ ~{threads} merged_~{fam_member}.bam ~{sep=' ' bams}

        echo "Sorting merged BAM"
        samtools sort -@ ~{threads} -o final_~{fam_member}.bam merged_~{fam_member}.bam

        echo "Indexing sorted BAM"
        samtools index -@ ~{threads} final_~{fam_member}.bam
    >>>

    output {
        File final = "final_~{fam_member}.bam"
        File final_idx = "final_~{fam_member}.bam.bai"
    }

    runtime {
        cpu: cpu
        gpu: false
        memory: "${mem}"
        disks: "local-disk ${disk_gb} SSD"
        docker: "${dockerSamtools}"
    }
}
