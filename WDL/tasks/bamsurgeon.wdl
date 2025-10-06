version 1.0

task bamsurgeon {
    input {
        File bam
        File bed
        File refGenomeBwaTar
        String fam_member
        String dockerBamsurgeon
    }

    # dynamic instance
    Int disk_gb = ceil( 2* (size(bam, "GiB"))) + 2
    String mem = "16 GB"
    Int threads = 16
    Int cpu = (threads)/2

    command <<<
        # debug
        echo "Running bamsurgeon with BAM: ~{bam} + BED: ~{bed} for the ~{fam_member}"

        tar -xvzf ~{refGenomeBwaTar}
        fasta=$(*.fasta)

        python3 /bamsurgeon/bin/addsnv.py \
            -v ~{bed} \
            -f ~{bam} \
            --aligner mem \
            --picardjar /picard.jar \
            -p 8 \
            -d 0.6 \
            -o ~{basename(bam)}.~{fam_member}.out.bam \
            -r $fasta

    >>>

    output {
        File spiked_bams = "~{basename(bam)}.~{fam_member}.out.bam"
    }

    runtime {
        cpu: cpu
        gpu: false
        memory: "${mem}"
        disks: "local-disk ${disk_gb} SSD"
        docker: "${dockerBamsurgeon}"
    }
}

task debug_print {
    input {
        File bam
        File bed
    }
    command <<<
        echo "BAM: ~{bam}"
        echo "BED: ~{bed}"
    >>>

    runtime {
        docker: "ubuntu:22.04"
    }
}