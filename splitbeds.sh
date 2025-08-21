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
    local progress=false
    local error=false
    local max_jobs=8
    local timeout=300  # force timeout after 5 mins 
    local running_containers=()

    for sample_dir in intermediates/*; do
        sample=$(basename "$sample_dir")
        beds=("$sample_dir"/*.bed)
        mkdir -p "outputs/$sample"

        for bedfile in "${beds[@]}"; do
            chrom=$(basename "$bedfile" .bed)
            outbam="outputs/$sample/${chrom}.bam"

            if [[ -f "$outbam" ]]; then
                echo "Skipping $outbam (already exists)"
                continue
            fi

            # launch container in detached mode and capture its container ID
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
            running_containers+=("$cid")
            progress=true

            # throttle: if max_jobs running, wait until at least one finishes/times out
            while (( ${#running_containers[@]} >= max_jobs )); do
                sleep 5
                # check status of each container
                for i in "${!running_containers[@]}"; do
                    cid="${running_containers[i]}"
                    status=$(docker inspect -f '{{.State.Status}}' "$cid")
                    if [[ "$status" == "exited" ]]; then
                        unset 'running_containers[i]'
                    fi
                done
            done
        done
    done

    # force kill jobs as they don't exit on their own
    for cid in "${running_containers[@]}"; do
        # wait up to 5 mins
        elapsed=0
        while [[ "$(docker inspect -f '{{.State.Status}}' "$cid")" != "exited" && "$elapsed" -lt "$timeout" ]]; do
            sleep 5
            elapsed=$((elapsed+5))
        done
        # kill
        if [[ "$(docker inspect -f '{{.State.Status}}' "$cid")" != "exited" ]]; then
            echo "Container $cid exceeded timeout, killing..."
            docker rm -f "$cid"
            error=true
        fi
    done

    return $error
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