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
    # track whether new bam is successfully created
    local progress=false
    # collect background container IDs to wait on
    local -a cids=()
    local error=false  

    # run bamsurgeon for each intermediate
    for sample_dir in intermediates/*; do     # father, mother, proband dirs
        sample=$(basename "$sample_dir")  
        beds=("$sample_dir"/*.bed)           # list all chromosome bed files for sample

        mkdir -p "outputs/$sample"

    for bedfile in "${beds[@]}"; do
            chrom=$(basename "$bedfile" .bed) # chromosomes (chr1, chr2.. etc)
            outbam="outputs/$sample/${chrom}.bam"

            # if chr{i}.bam exists then skip 
            if [[ -f "$outbam" ]]; then
                echo "Skipping $outbam (already exists)"
                continue
            fi

            # launch detached, record container id, and mark progress
            progress=true
            cid=$(docker run --rm -d -v "$(pwd)":/data bamsurgeon-env \
                python3 /bamsurgeon/bin/addsnv.py \
                    -v "/data/$bedfile" \
                    -f "/data/bam/WGS_EX2500218_22CFV7LT4.bam" \
                    --aligner mem \
                    --picardjar /picard.jar \
                    -p 8 \
                    -d 0.6 \
                    -o "/data/outputs/$sample/${chrom}.bam" \
                    -r "$REFERENCE_GENOME" </dev/null)
            cids+=("$cid")
         done       
    done

    # wait for all launched containers to complete and capture failures
    if [ ${#cids[@]} -gt 0 ]; then
        for cid in "${cids[@]}"; do
            status=$(docker wait "$cid" || true)
            if [ "${status:-1}" != "0" ]; then
                error=true
                echo "Container $cid exited with status $status" >&2
                # best-effort logs (container may be auto-removed on exit with --rm)
                docker logs "$cid" 2>&1 | sed "s/^/[${cid}] /" >&2 || true
            fi
        done
    fi
    if [ "$progress" = true ] && [ "$error" = false ]; then
        return 0
    fi
    return 1
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