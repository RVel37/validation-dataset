version 1.0

import "tasks/bam_to_fastq.wdl" as bam_to_fastq

struct SampleInputs {
    File bam
    String fam_member
    String proband_sex
}

workflow fastqs {
    input {
        Array[SampleInputs] samples
        String dockerHtslib
    }

    scatter (s in samples) {
        call bam_to_fastq.bam_to_fastq {
            input:
            bam = s.bam,
            fam_member=s.fam_member,
            dockerHtslib=dockerHtslib
        }

        call bam_to_fastq.name{
            r1_fastq = bam_to_fastq.r1_fastq
            r2_fastq = bam_to_fastq.r2_fastq
            sex = proband_sex
        }
    }

    output {
        # fastqs
        Array[File] r1_fastqs = bam_to_fastq.r1_fastq
        Array[File] r2_fastqs = bam_to_fastq.r2_fastq
    }
}