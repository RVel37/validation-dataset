version 1.0

task bam_to_fastq {
    input {
        File bam
        String fam_member
        String dockerHtslib
    }

    Int disk_gb = ceil(3 * size(bam, "GiB"))
    String mem = "16 GB"
    Int threads = 16
    Int task_threads = 8
    Int cpu = (threads)/2

    command <<<

    # convert coord -> queryname sorted bam
    samtools sort -n \
    -@ ~{task_threads} \
    -o qsorted_~{fam_member}.bam \
    ~{bam}

    # fixmate (rebuild mate fields & correct flags)
    samtools fixmate \
    -@ ~{task_threads} \
    -m qsorted_~{fam_member}.bam \
    fixmate_~{fam_member}.bam

    # convert to bgzipped fastq (files automatically compressed)
    samtools fastq -@ 8 \
        -1 R1_~{fam_member}.fastq.gz \
        -2 R2_~{fam_member}.fastq.gz \
        fixmate_~{fam_member}.bam

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

task name {
    input {
        File r1_fastq
        File r2_fastq
        File proband_sex
    }

    Int disk_gb = ceil(1.2 * size((r1_fastq) + (r2_fastq), "GiB"))
    String mem = "4 GB"
    Int threads = 4
    Int cpu = (threads)/2

    command <<<
    
    # Determine sample no. based on filename
    if [[ "~{basename(r1_fastq)}" == *proband* ]]; then
        s_no="~{S1}"
    elif [[ "~{basename(r1_fastq)}" == *mother* ]]; then
        s_no="~{S2}"
    elif [[ "~{basename(r1_fastq)}" == *father* ]]; then
        s_no="~{S3}"
    fi

    # Build final output names
    R1_name="WGS_EX_~{proband_sex}_123_${s_no}_L001_R1_001.fastq.gz"
    R2_name="WGS_EX_~{proband_sex}_123_${s_no}_L001_R2_001.fastq.gz"

    mv "~{r1_fastq}" "${R1_name}"
    mv "~{r2_fastq}" "${R2_name}"

    >>>

    output {
        File namedR1 = ${R1_name}
        File namedR2 = ${R2_name}
    }

    runtime {
        cpu: cpu
        gpu: false
        memory: "${mem}"
        disks: "local-disk ${disk_gb} SSD"
        docker: "${dockerHtslib}"
    }

}