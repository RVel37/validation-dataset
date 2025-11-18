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
    String mem = "16 GB"
    Int threads = 8
    Int cpu = (threads)/2

    command <<<
        echo "usage at start ($(date))"; free -h; df -h /
        du -h "~{bam}"

        echo "Running bamsurgeon with BAM: ~{bam} + BED: ~{bed} for the ~{fam_member}"

        # unpack reference genome
        mkdir -p ref
        tar -zxvf ~{refGenomeBwaTar} -C ref --no-same-owner
        referenceFasta=$(ls ref/*.fasta | head -n1)

        outbam="~{basename(bam, ".bam")}.~{fam_member}.out.bam"

        timeout 3600 \ # bamsurgeon occasionally gets stuck if no mutations are added
        python3 /usr/local/bin/addsnv.py \
                -v ~{bed} \
                -f ~{bam} \
                -r ${referenceFasta} \
                --aligner mem \
                --picardjar /usr/local/bin/picard.jar \
                -p 8 \
                -d 0.2 \
                --force \
                -o ${outbam}

        echo "usage at end ($(date))"; free -h; df -h /

        if [ ! -s "${outbam}" ] ; then
            echo "No spiked bam was generated for this chromosome. Using original bam. "
            # if addsnv produced nothing, copy input bam instead
            cp "~{bam}" "${outbam}"
        fi

    >>>

    output {
        File spiked_bams = "~{basename(bam, '.bam')}.~{fam_member}.out.bam"
    }

    runtime {
        cpu: cpu
        gpu: false
        memory: "${mem}"
        disks: "local-disk ${disk_gb} SSD"
        docker: "${dockerBamsurgeon}"
    }
}
