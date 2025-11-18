version 1.0

import "tasks/split_samples.wdl" as split_samples
import "tasks/bamsurgeon.wdl" as bamsurgeon
import "tasks/merge_bams.wdl" as merge_bams
import "tasks/bam_to_fastq.wdl" as bam_to_fastq

struct SampleInputs {
    File bam
    File bai
    File bed
    String fam_member
}

workflow main {
    input {
        Array[SampleInputs] samples
        File refGenomeBwaTar
        String dockerSamtools
        String dockerBamsurgeon
        String dockerHtslib
    }

    scatter (s in samples) {
        call split_samples.split_samples as split {
            input:
                bam = s.bam,
                bai = s.bai,
                bed = s.bed,
                fam_member=s.fam_member,
                dockerSamtools=dockerSamtools
        }

        # scatter over aligned BED and BAM arrays
        scatter (i in range(length(split.bam_array))) {

            call bamsurgeon.bamsurgeon as spike_in {
                input:
                    bam = split.bam_array[i],
                    bai = split.bai_array[i],
                    bed = split.bed_array[i],
                    refGenomeBwaTar = refGenomeBwaTar,
                    fam_member = s.fam_member,
                    dockerBamsurgeon = dockerBamsurgeon
            }
        }

        # collect the output of bamsurgeon task
        Array[File] spiked_bam_array = spike_in.spiked_bams

        call merge_bams.merge_bams {
            input:
            bams = spiked_bam_array,
            fam_member=s.fam_member,
            dockerSamtools=dockerSamtools
        }

        call bam_to_fastq.bam_to_fastq {
            input:
            bam = merge_bams.bam,
            fam_member=s.fam_member,
            dockerHtslib=dockerHtslib
        }

    }

    output {
        # collect coord-sorted bams (enables inspection in IGV)
        Array[File] final_bams = merge_bams.bam
        Array[File] final_bais = merge_bams.bai

        # fastqs (bam to fastq step, can be executed in separate fastq workflow)
        Array[File] r1_fastqs = bam_to_fastq.r1_fastq
        Array[File] r2_fastqs = bam_to_fastq.r2_fastq
    }
}