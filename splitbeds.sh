#!/bin/bash
set -euo pipefail
 
######### FUNCTION DEFINITIONS ###########
 
 
# split each bed by chromosome
split_beds() {
    for bedfile in ./*.bed; do
        sample=$(basename "$bedfile" .bed)   # father, mother, proband
        mkdir -p "intermediates/$sample"
 
        awk -v outdir="intermediates/$sample" \
            '{ print > (outdir "/" $1 ".bed") }' "$bedfile"
    done
}
 
# create small bam files for each intermediate bed
create_bams() {
    # variable to prompt function to keep going (otherwise stops after creating a single bam file)
    local progress=false  
 
    # run bamsurgeon for each intermediate
    for i in intermediates/*/*.bed; do
        sample=$(basename "$(dirname "$i")")  # father, mother, proband
        chrom=$(basename "$i" .bed)           # chromosomes (chr1, chr2.. etc)
        outbam="outputs/$sample/${chrom}.bam" # define output bam
 
        mkdir -p "outputs/$sample"
 
        # if num.bam exists then skip creation of that one
        if [[ -f "$outbam" ]]; then
            echo "Skipping $outbam (already exists)"
            continue
        fi
        docker run --rm -d -v "$(pwd)":/data bamsurgeon-env \
            python3 /bamsurgeon/bin/addsnv.py \
                -v "/data/$i" \
                -f "/data/bam/WGS_EX2500218_22CFV7LT4.bam" \
                --aligner mem \
                --picardjar /picard.jar \
                -p 8 \
                -d 0.6 \
                -o "/data/outputs/$sample/${chrom}.bam" \
                -r "/data/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_hs38d1_maskedGRC_exclusions_v2_no_chr.fasta" </dev/null
        progress=true    
    done
    $progress && return 0 || return 1
}
 
# merge bam files with samtools
merge_bams() {
    for sample in father mother proband; do
        echo "Merging BAMs for $sample..."
        docker run --rm -v "$(pwd)":/data bamsurgeon-env \
            samtools merge -@ 8 "/data/outputs/${sample}/${sample}_merged.bam" \
            /data/outputs/${sample}/*.bam
 
        docker run --rm -v "$(pwd)":/data bamsurgeon-env \
            samtools index "/data/outputs/${sample}/${sample}_merged.bam"
    done
}
 
 
######### MAIN EXECUTION FLOW ###########
split_beds
 
# keep running create_bams until all are finished
while create_bams; do
    echo "Checking for remaining bed files without matching bams..."
done
 
# verify that all the bams have beds
for i in father mother proband; do
    bed_count=$(find "intermediates/$i" -type f -name "*.bed" | wc -l)
    # don't count any "_merged.bam" files we may have previously created
    bam_count=$(find "outputs/$i" -type f -name "*.bam" ! -name "*_merged.bam" | wc -l)
 
    if [[ "$bed_count" -ne "$bam_count" ]]; then
        echo "ERROR: $i is missing some BAM files." >&2
        exit 1
    else
        echo "All intermediate bams for the $i sample created. Merging..."
    fi
done
 
# now merge everything
merge_bams