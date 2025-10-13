version 1.0

task merge_bams {
    input {
        Array[File?] bams
        String fam_member
        String dockerSamtools
    }

    Int disk_gb = 64
    String mem = "16 GB"
    Int threads = 8
    Int cpu = 8

    Array[File] present_bams = select_all(bams) # remove nulls

    command <<<

        echo "List of BAMs:"
        for bam in ~{sep=' ' present_bams}; do
            echo $bam
        done

        echo "Merging BAMs for ~{fam_member}"
        samtools merge -@ 8 merged_~{fam_member}.bam ~{sep=' ' present_bams}

        echo "Sorting merged BAM"
        samtools sort -@ 8 -o final_~{fam_member}.bam merged_~{fam_member}.bam

        echo "Indexing sorted BAM"
        samtools index -@ 8 final_~{fam_member}.bam
    >>>

    output {
        File? final = "final_~{fam_member}.bam"
        File? final_idx = "final_~{fam_member}.bam.bai"
    }

    runtime {
        cpu: cpu
        gpu: false
        memory: "${mem}"
        disks: "local-disk ${disk_gb} SSD"
        docker: "${dockerSamtools}"
    }
}
