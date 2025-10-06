version 1.0

import "tasks/split_samples.wdl" as split_samples
import "tasks/bamsurgeon.wdl" as bamsurgeon
import "tasks/merge_bams.wdl" as merge_bams

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
            call debug_print {
                input:
                    bam = split.bam_array[i],
                    bed = split.bed_array[i],
            }

            call bamsurgeon.bamsurgeon as spike {
                input:
                    bam = split.bam_array[i],
                    bed = split.bed_array[i],
                    refGenomeBwaTar = refGenomeBwaTar,
                    fam_member = s.fam_member,
                    dockerBamsurgeon = dockerBamsurgeon
            }
        }

        # collect the output of bamsurgeon task
        Array[File] spiked_bam_array = spike.spiked_bams

        call merge_bams.merge_bams {
            input:
            bams = spiked_bam_array,
            fam_member=s.fam_member
            dockerSamtools=dockerSamtools
        }
    }

    output {
        Array[File] bams = merge_bams.final
        Array[File] bais = merge_bams.final_idx
    }
}