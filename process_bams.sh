#!/bin/bash
set -euo pipefail
mkdir -p outputs split_beds split_bams

REFERENCE_GENOME="/data/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_hs38d1_maskedGRC_exclusions_v2_no_chr.fasta"

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

    # if [[ ! -f "$bamfile.bai" ]]; then
    #     samtools index "$bamfile"
    # fi

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
    max_jobs=3

    # collect all beds that still need a bam
    beds_to_process=()
    for sample_dir in intermediates/*; do
        sample=$(basename "$sample_dir")

        for bedfile in "$sample_dir"/*.bed; do
            chrom=$(basename "$bedfile" .bed)
            outbam="outputs/$sample/${chrom}.bam"
            
            if [[ ! -f "$outbam" ]]; then
                beds_to_process+=("$bedfile")
            fi
        done
    done

    # if nothing left, exit
    if [[ ${#beds_to_process[@]} -eq 0 ]]; then
    echo "No beds to process."
        return 1  
    fi

    # process in batches
    i=0
    # while the list of jobs is less than total number of items in the list (beds)
    while [[ $i -lt ${#beds_to_process[@]} ]]; do
        # take up to 3 items and process this as a batch
        batch=("${beds_to_process[@]:i:max_jobs}")
        running_containers=()

        # launch batch
        for bedfile in "${batch[@]}"; do
            sample=$(basename "$(dirname "$bedfile")")
            chrom=$(basename "$bedfile" .bed)
            outbam="outputs/$sample/${chrom}.bam"

            cid=$(docker run -d -v "$(pwd)":/data bamsurgeon-env \
                python3 /bamsurgeon/bin/addsnv.py \
                    -v "/data/$bedfile" \
                    -f "/data/split_bams/${sample}/${chrom}.bam" \
                    --aligner mem \
                    --picardjar /picard.jar \
                    -p 8 \
                    -d 0.6 \
                    -o "/data/outputs/$sample/${chrom}.bam" \
                    -r "$REFERENCE_GENOME") 2>&1
            echo "$(date '+%H:%M:%S') Started $cid for $sample sample, chr $chrom"
            running_containers+=("$cid")
        done
            wait # for all the jobs in the batch to finish
            i=$((i + max_jobs))
        done
    return 0
}


# merge bam files with samtools
merge_bams() {
    for sample_dir in outputs/*; do
        samples=$(basename "$sample_dir")
        for sample in samples; do 
        docker run -v "$(pwd)":/data biocontainers/samtools:v1.9-4-deb_cv1 \
            samtools merge $sample.full.bam /data/$sample/* 

}


### EXECUTION WORKFLOW ###

split_beds
split_bams
spike_in_variants

# verify that all the bams have beds 
#   (this is currently flagging that no mitochondrial variants are being processed)
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
#merge_bams
