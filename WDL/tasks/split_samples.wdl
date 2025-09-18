version 1.0


task split_samples {
    input {
        File bam
        File bai
        File bed
        String fam_member
    }

    # dynamic instance
    Int disk_gb = ceil( 2* (size(bam, "GiB"))) + 2
    String mem = "8 GB"
    Int threads = 8
    Int cpu = (threads)/2

    command <<<

        sample=~{fam_member}


        # Split BED by chromosome
        awk -v outdir="split_beds/$sample" '{ print > (outdir "/" $1 ".bed") }' "~{bed}"

        # Split BAM by chromosome + index
        for chrom in $(samtools idxstats "~{bam}" | cut -f1 | grep -v '\*' | sort -V); do
            bam_out="split_bams/$sample/${chrom}.bam"
            bed_out="split_beds/$sample/${chrom}.bed"
            samtools view -b "~{bam}" "$chrom" > "$bam_out"
            samtools index "$bam_out"
        done

    >>>

    output {
        Array[File] bam_array = sort(glob("split_bams/~{fam_member}/*.bam"))
        Array[File] bed_array = sort(glob("split_beds/~{fam_member}/*.bed"))
    }

    runtime {

    }
}


struct Chromosome {
    File chr_bam
    File chr_bed
    String chrom
}

task pair_chromosomes {
    input {
        Array[File] bam_array
        Array[File] bed_array
    }
    
    output {
        Array[Chromosome] Chr = [
            for (i in range(length(bam_array))) {
                Chromosome(
                    chr_bam = bam_array[i],
                    chr_bed = bed_array[i],
                    chrom   = basename(bam_array[i], ".bam")
                )
            }
        ]
    }
}
