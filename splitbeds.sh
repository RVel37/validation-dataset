#!/bin/bash
set -euo pipefail
mkdir -p outputs intermediates

REFERENCE_GENOME="/data/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_hs38d1_maskedGRC_exclusions_v2_no_chr.fasta"

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
    # local variables - won't get overwritten when running tasks concurrently
    local max_jobs=3
    local timeout=180 # 3 mins
    local progress=false

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
                    -f "/data/bams/${sample}.bam" \
                    --aligner mem \
                    --picardjar /picard.jar \
                    -p 8 \
                    -d 0.6 \
                    -o "/data/outputs/$sample/${chrom}.bam" \
                    -r "$REFERENCE_GENOME")
            echo "$(date '+%H:%M:%S') Started $cid for $sample sample, chr $chrom"
            running_containers+=("$cid")
        
        # count how long individual task runs for; terminate after 3 mins
        (
            elapsed=0
            while [[ "$(docker inspect -f '{{.State.Status}}' "$cid")" != "exited" && "$elapsed" -lt "$timeout" ]]; do
                sleep 5
                elapsed=$((elapsed+5))
            done
            if [[ "$(docker inspect -f '{{.State.Status}}' "$cid")" != "exited" ]]; then
                echo "$(date '+%H:%M:%S') $cid has been running for 3 mins, killing task. "
                docker rm -f "$cid"
            fi
        ) &
    done
        wait # for all the jobs in the batch to finish
        i=$((i + max_jobs))
    done

    return 0
}


# merge bam files with samtools
merge_bams() {
    for sample in father mother proband; do
        echo "Merging BAMs for $sample..."
        docker run --rm -d -v "$(pwd)":/data bamsurgeon-env \
            samtools merge -@ 8 "/data/outputs/${sample}/${sample}_merged.bam" \
            /data/outputs/${sample}/*.bam

        docker run --rm -v "$(pwd)":/data bamsurgeon-env \
            samtools index "/data/outputs/${sample}/${sample}_merged.bam"
    done
}


######### MAIN EXECUTION FLOW ###########
split_beds

create_bams

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