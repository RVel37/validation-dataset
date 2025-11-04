version 1.0

task bam_to_fastq {
    input {
        File coord_bam
        String fam_member
        String dockerSamtools
    }

    Int disk_gb = ceil(1.5 * size(coord_bam, "GiB"))
    String mem = "16 GB"
    Int threads = 16
    Int cpu = (threads)/2
	Int task_threads = (threads) - 4

    command <<<
	# fixmate (rebuild mate fields and correct mate flags)
	samtools fixmate \
	-@ ~{task_threads} \
	-m ~{coord_bam} \
	fixmate_~{fam_member}.bam

	# bam to fastq
	samtools fastq -@ ~{task_threads} \
	-1 R1_~{fam_member}.fastq \
	-2 R2_~{fam_member}.fastq \
	-0 /dev/null \
	-s /dev/null \
	fixmate_~{fam_member}.bam

	>>>

    output {
        File r1_fastq = "R1_~{fam_member}.fastq"
        File r2_fastq = "R2_~{fam_member}.fastq"
    }

    runtime {
        cpu: cpu
        gpu: false
        memory: "${mem}"
        disks: "local-disk ${disk_gb} SSD"
        docker: "${dockerSamtools}"
    }
}
