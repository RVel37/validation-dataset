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
    Int disk_gb = ceil( 2* (size(bam, "GiB"))) + 2
    String mem = "32 GB"
    Int threads = 16
    Int cpu = (threads)/2

    command <<<
        
        echo "Running bamsurgeon with BAM: ~{bam} + BED: ~{bed} for the ~{fam_member}"

        # unpack reference genome
        mkdir -p ref
        tar -zxvf ~{refGenomeBwaTar} -C ref --no-same-owner
        referenceFasta=$(ls ref/*.fasta | head -n1)
        
        echo "DEBUG: Checking for picard.jar"
        find / -maxdepth 5 -type f -name "picard.jar" 2>/dev/null || echo "picard.jar not found!"


        python3 /usr/local/bin/addsnv.py \
                -v ~{bed} \
                -f ~{bam} \
                -r ${referenceFasta} \
                --aligner mem \
                --picardjar /usr/local/bin/picard.jar \
                -p 12 \
                -d 0.6 \
                -o ~{basename(bam)}.~{fam_member}.out.bam 

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


#########################################################################

task debug_print {
    input {
        File bam
        File bai
        File bed
    }
    command <<<
        echo "BAM: ~{bam}"
        echo "BAI: ~{bai}"
        echo "BED: ~{bed}"
    >>>

    runtime {
        docker: "ubuntu:24.04"
    }
}