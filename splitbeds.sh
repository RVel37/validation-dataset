#!/bin/bash
# Split each bed by chromosome
for bedfile in ./*.bed; do
    sample=$(basename "$bedfile" .bed)   # father, mother, proband
    mkdir -p "intermediates/$sample"

    awk -v outdir="intermediates/$sample" \
        '{ print > (outdir "/" $1 ".bed") }' "$bedfile"
done

# Run bamsurgeon for each intermediate
for i in intermediates/*/*.bed; do
    sample=$(basename "$(dirname "$i")")  # father, mother, proband
    chrom=$(basename "$i" .bed)           # chromosomes (chr1, chr2.. etc)

    mkdir -p "outputs/$sample"

    # if num.bam exists then skip creation of that one

    docker run --rm -v "$(pwd)":/data bamsurgeon-env \
        python3 /bamsurgeon/bin/addsnv.py \
            -v "/data/$i" \
            -f "/data/bam/WGS_EX2500218_22CFV7LT4.bam" \
            --aligner mem \
            --picardjar /picard.jar \
            -p 8 \
            -o "/data/outputs/$sample/${chrom}.bam" \
            -r "/data/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_hs38d1_maskedGRC_exclusions_v2_no_chr.fasta"
done

# # merge bam files with samtools
# for sample in father mother proband; do

#     docker run --rm -v "$(pwd)":/data bamsurgeon-env \
#         samtools merge -@ 8 "/data/outputs/${sample}/${sample}_merged.bam" \
#         /data/outputs/${sample}/*.bam
#     docker run --rm -v "$(pwd)":/data bamsurgeon-env \
#         samtools index "/data/outputs/${sample}/${sample}_merged.bam"
# done