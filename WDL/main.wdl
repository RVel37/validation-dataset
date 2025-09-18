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


        scatter (i in range(length(split.bam_array))) {
            call bamsurgeon.bamsurgeon {
                input:
                    bam = split.bam_array[i],
                    bed = split.bed_array[i],
                    refGenomeBwaTar = refGenomeBwaTar,
                    fam_member = s.fam_member,
                    dockerBamsurgeon = dockerBamsurgeon
            }
        }

        bams = flatten(bamsurgeon.spiked_bams)

        call merge_bams.merge_bams {
            input:
            bams = bamsurgeon.spiked_bams,
            fam_member=s.fam_member
        }
    }

    output {
        Array[File] bams = merge_bams.merged
    }
}