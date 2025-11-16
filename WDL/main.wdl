version 1.0

import "tasks/bam_to_fastq.wdl" as bam_to_fastq

struct SampleInputs {
    File bam
    String fam_member
}

workflow main {
    input {
        Array[SampleInputs] samples
        String dockerHtslib
    }

    scatter (s in samples) {
        call bam_to_fastq.bam_to_fastq {
            input:
            coord_bam = s.bam,
            fam_member=s.fam_member,
            dockerHtslib=dockerHtslib
            
        }    
    }

    output {
        # fastqs
        Array[File] r1_fastqs = bam_to_fastq.r1_fastq
        Array[File] r2_fastqs = bam_to_fastq.r2_fastq
    }
}