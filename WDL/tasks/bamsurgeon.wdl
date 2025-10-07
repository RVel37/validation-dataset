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
        
        echo "Running bamsurgeon with BAM: ~{bam} + BED: ~{bed} for the ~{fam_member}"

        # unpack reference genome
        mkdir -p ref
        tar -zxvf ~{refGenomeBwaTar} -C ref --no-same-owner
        referenceFasta=$(ls ref/*.fasta | head -n1)
        
        echo "ls:"; ls; ls ref

        # DEBUG - CANT FIND PATH #
        ADDSNV_PATH=$(find /bamsurgeon /usr /opt / -type f -name addsnv.py 2>/dev/null | head -n 1)

        if [[ -z "$ADDSNV_PATH" ]]; then
            echo "ERROR: addsnv.py not found anywhere on system!" >&2
            echo "Contents of /bamsurgeon/bin (if exists):"
            ls -l /bamsurgeon/bin 2>/dev/null || echo "No /bamsurgeon/bin directory."
            exit 1
        fi

        echo "Found addsnv.py at: $ADDSNV_PATH"

        python3 "$ADDSNV_PATH" \
                -v ~{bed} \
                -f ~{bam} \
                -r ${referenceFasta} \
                --aligner mem \
                --picardjar /picard.jar \
                -p 8 \
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
        docker: "ubuntu:24.04"
    }
}