#!/bin/bash
set -euo pipefail
mkdir -p outputs split_beds split_bams

REFERENCE_GENOME="/reference_genome.fasta"

# Pull required Docker images
docker pull biocontainers/samtools:v1.9-4-deb_cv1
docker build -t bamsurgeon-env .

# split each bed by chromosome
split_beds() {
    for bedfile in ./*.bed; do
        # create basenames father, mother, proband
        sample=$(basename "$bedfile" .bed)   
        mkdir -p "split_beds/$sample"

        # create output directories for each
        awk -v outdir="split_beds/$sample" \
            '{ print > (outdir "/" $1 ".bed") }' "$bedfile"
    done
    # split_beds/family member/*.bed
}


# split each bam by chromosome
split_bams() {
    mkdir -p split_bams
    for bamfile in bams/*.bam; do
        sample=$(basename "$bamfile" .bam)
        mkdir -p "split_bams/$sample"

    for chrom in $(samtools idxstats "$bamfile" | cut -f1 | grep -v '\*'); do
        outbam="split_bams/$sample/${chrom}.bam"
        if [[ ! -f "$outbam" ]]; then
            echo "Splitting $sample $chrom"
            samtools view -b "$bamfile" "$chrom" > "$outbam"
            samtools index "$outbam"
            # split_bams/family member/*.bam
        fi
    done
done
}


# run bamsurgeon for all intermediates
spike_in_variants() {
    for bedfile in split_beds/*/*.bed; do
        sample=$(basename "$(dirname "$bedfile")")
        chrom=$(basename "$bedfile" .bed)
        mkdir -p "outputs/$sample"
        outbam="outputs/$sample/${chrom}.bam"

        echo "$(date '+%H:%M:%S') Running bamsurgeon on $sample $chrom"
        docker run --rm -v "$(pwd)":/data bamsurgeon-env \
            python3 /bamsurgeon/bin/addsnv.py \
                -v "/data/$bedfile" \
                -f "/data/split_bams/${sample}/${chrom}.bam" \
                --aligner mem \
                --picardjar /picard.jar \
                -p 8 \
                -d 0.6 \
                -o "/data/$outbam" \
                -r "$REFERENCE_GENOME"
    done
}


# merge bam files with samtools after - dont do this yet
merge_bams() {
    for sample_dir in outputs/*; do
        sample=$(basename "$sample_dir")
        echo "Merging BAMs for $sample"
        
        docker run --rm -v "$(pwd)":/data biocontainers/samtools:v1.9-4-deb_cv1 \
            samtools merge \
            "/data/outputs/${sample}.full.bam" \
            /data/outputs/"$sample"/*.bam
    done
}


### EXECUTION WORKFLOW ###

split_beds
split_bams
spike_in_variants

# verify that all the bams have beds 
#   (BUG: this is currently flagging that no mitochondrial variants are being processed)
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
