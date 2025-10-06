version 1.0

task merge_bams {
    input {
        Array[File] bams
        String fam_member
        String dockerSamtools
    }

    Int disk_gb = ceil(2 * size(bams, "GiB")) + 2
    String mem = "16 GB"
    Int threads = 16
    Int cpu = (threads)/2

    command <<<
        echo "Merging BAMs for ~{fam_member}:"
        samtools merge -@ ~{threads} merged_~{fam_member}.bam ~{sep=' ' bams}

        echo "Sorting merged BAM:"
        samtools sort -@ ~{threads} -o final_~{fam_member}.bam merged_~{fam_member}.bam

        echo "Indexing sorted BAM:"
        samtools index final_~{fam_member}.bam
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
