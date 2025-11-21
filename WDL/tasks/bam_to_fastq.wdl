version 1.0

task bam_to_fastq {
    input {
        File bam
        String fam_member
        String dockerSamtools
    }

    Int disk_gb = ceil(4 * size(bam, "GiB"))
    String mem = "16 GB"
    Int threads = 16
    Int task_threads = 12
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
        docker: "${dockerSamtools}"
    }
}


task name {
    input {
        File r1_fastq
        File r2_fastq
        String fam_member
        String proband_sex
    }

    Int disk_gb = ceil(1.2 * ( size(r1_fastq, "GiB") + size(r2_fastq, "GiB") ))
    String mem = "4 GB"
    Int threads = 4
    Int cpu = (threads)/2

    command <<<
    
    # Determine sample no. based on filename
    if [[ "~{fam_member}" == "proband" ]]; then
        s_no="S1"
    elif [[ "~{fam_member}" == "mother" ]]; then
        s_no="S2"
    elif [[ "~{fam_member}" == "father" ]]; then
        s_no="S3"
    fi

    # Build final output names
    R1_name="WGS_EX~{fam_member}~{proband_sex}_123_${s_no}_L001_R1_001.fastq.gz"
    R2_name="WGS_EX~{fam_member}~{proband_sex}_123_${s_no}_L001_R2_001.fastq.gz"

    mv ~{r1_fastq} "${R1_name}"
    mv ~{r2_fastq} "${R2_name}"
    >>>

    output {
        File namedR1 = glob("WGS_EX~{fam_member}~{proband_sex}_123_S*_L001_R1_001.fastq.gz")[0]
        File namedR2 = glob("WGS_EX~{fam_member}~{proband_sex}_123_S*_L001_R2_001.fastq.gz")[0]
    }

    runtime {
        cpu: cpu
        gpu: false
        memory: "${mem}"
        disks: "local-disk ${disk_gb} SSD"
        docker: "ubuntu:latest"
    }
}
