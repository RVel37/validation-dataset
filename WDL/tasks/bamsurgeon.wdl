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
    String mem = "8 GB"
    Int threads = 4
    Int cpu = (threads)/2

    command <<<

        tar -xvzf ~{refGenomeBwaTar}
        fasta=$(*.fasta)

        python3 /bamsurgeon/bin/addsnv.py \
            -v ~{bed} \
            -f ~{bam} \
            --aligner mem \
            --picardjar /picard.jar \
            -p 8 \
            -d 0.6 \
            -o ~{basename(bam)}.{fam_member}.out.bam \
            -r $fasta

    >>>

    output {
        File spiked_bams = "~{basename(bam)}.{fam_member}.out.bam"
    }

    runtime {
        cpu: cpu
        gpu: false
        memory: "${mem}"
        disks: "local-disk ${disk_gb} SSD"
        docker: {dockerBamsurgeon}
    }
}