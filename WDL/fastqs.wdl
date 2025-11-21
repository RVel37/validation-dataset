version 1.0

import "tasks/bam_to_fastq.wdl" as bam_to_fastq

struct SampleInputs {
    File bam
    String fam_member
}

workflow fastq {
    input {
        Array[SampleInputs] samples
        String dockerSamtools
        String proband_sex
    }

    scatter (s in samples) {
        call bam_to_fastq.bam_to_fastq {
            input:
            bam = s.bam,
            fam_member=s.fam_member,
            dockerSamtools=dockerSamtools
        }

        call bam_to_fastq.name {
            input:
            r1_fastq = bam_to_fastq.r1_fastq,
            r2_fastq = bam_to_fastq.r2_fastq,
            fam_member=s.fam_member,
            proband_sex = proband_sex
        }
    }

    output {
        Array[File] r1_fastqs = bam_to_fastq.r1_fastq
        Array[File] r2_fastqs = bam_to_fastq.r2_fastq
        # fastqs
        Array[File] r1_named = name.namedR1
        Array[File] r2_named = name.namedR2
    }
}
