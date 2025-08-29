#!/bin/bash

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
                    -r "$REFERENCE_GENOME")>file.log 2>&1
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

create_bams