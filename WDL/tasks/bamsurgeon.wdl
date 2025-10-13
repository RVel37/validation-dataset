version 1.0

task bamsurgeon {
    input {
        File bam
        File bai
        File bed
        File refGenomeBwaTar
        String fam_member
        String dockerBamsurgeon
    }

    # dynamic instance
    Int disk_gb = ceil((size(bam, "GiB"))) + 10
    String mem = "32 GB"
    Int threads = 8
    Int cpu = (threads)

    command <<<
        echo "usage at start ($(date))"; free -h; df -h /
        du -h "~{bam}"

        echo "Running bamsurgeon with BAM: ~{bam} + BED: ~{bed} for the ~{fam_member}"

        # unpack reference genome
        mkdir -p ref
        tar -zxvf ~{refGenomeBwaTar} -C ref --no-same-owner
        referenceFasta=$(ls ref/*.fasta | head -n1)

        python3 /usr/local/bin/addsnv.py \
                -v ~{bed} \
                -f ~{bam} \
                -r ${referenceFasta} \
                --aligner mem \
                --picardjar /usr/local/bin/picard.jar \
                -p 8 \
                -d 0.6 \
                -o ~{basename(bam, ".bam")}.~{fam_member}.out.bam

        echo "usage at end ($(date))"; free -h; df -h /

    >>>

    output {
        File? spiked_bams = "~{basename(bam, '.bam')}.~{fam_member}.out.bam"
    }

    runtime {
        cpu: cpu
        gpu: false
        memory: "${mem}"
        disks: "local-disk ${disk_gb} SSD"
        docker: "${dockerBamsurgeon}"
    }
}
