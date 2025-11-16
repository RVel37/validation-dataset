version 1.0

task bam_to_fastq {
    input {
        File coord_bam
        String fam_member
        String dockerHtslib
    }

    Int disk_gb = ceil(3 * size(coord_bam, "GiB"))
    String mem = "16 GB"
    Int threads = 16
    Int cpu = (threads)/2

    command <<<
	# bam to fastq

        samtools fastq -@ 8 \
        -1 >(bgzip -@ 4 > R1_~{fam_member}.fastq.gz) \
        -2 >(bgzip -@ 4 > R2_~{fam_member}.fastq.gz) \
        ~{coord_bam}
    
	>>>

    output {
        File r1_fastq = "R1_~{fam_member}.fastq.gz"
        File r2_fastq = "R2_~{fam_member}.fastq.gz"
    }

    runtime {
        cpu: cpu
        gpu: false
        memory: "${mem}"
        disks: "local-disk ${disk_gb} SSD"
        docker: "${dockerHtslib}"
    }
}
