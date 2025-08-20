#!/bin/bash
set -euo pipefail
mkdir -p outputs intermediates

######### FUNCTION DEFINITIONS ###########
 
# split each bed by chromosome
split_beds() {
    for bedfile in ./*.bed; do
        # create basenames father, mother, proband
        sample=$(basename "$bedfile" .bed)   
        mkdir -p "intermediates/$sample"

        # create output directories for each
        awk -v outdir="intermediates/$sample" \
            '{ print > (outdir "/" $1 ".bed") }' "$bedfile"
    done
}


# create bams corresponding to each intermediate bed
create_bams() {
    # track whether new bam is successfully created
    progress=false  

    # Loop over each sample (father, mother, proband)
    for sample_dir in intermediates/*; do
        sample=$(basename "$sample_dir")
        mkdir -p "outputs/$sample"

        # Get all chromosome BED files for this sample
        beds=( "$sample_dir"/*.bed )
        total=${#beds[@]} # total number of intermediate beds for sample

        batch_size=4
        index=0

        # Process BED files in batches
        while [ $index -lt $total ]; do
            # Take a batch of up to batch_size files
            batch=( "${beds[@]:$index:$batch_size}" )

            # Process each BED in the batch sequentially
            for bedfile in "${batch[@]}"; do
                chrom=$(basename "$bedfile" .bed)
                outbam="outputs/$sample/${chrom}.bam"

                if [[ -f "$outbam" ]]; then
                    echo "Skipping $outbam (already exists)"
                    continue
                fi
                
                docker run --rm -v "$(pwd)":/data bamsurgeon-env \
                    python3 /bamsurgeon/bin/addsnv.py \
                        -v "/data/$bedfile" \
                        -f "/data/bams/${sample}.bam" \
                        --aligner mem \
                        --picardjar /picard.jar \
                        -p 8 \
                        -d 0.6 \
                        -o "/data/outputs/$sample/${chrom}.bam" \
                        -r "/data/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_hs38d1_maskedGRC_exclusions_v2_no_chr.fasta" </dev/null

                progress=true
            done

            # Move to the next batch
            index=$((index + batch_size))
        done
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

# keep running create_bams until all are finished - retry if not
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