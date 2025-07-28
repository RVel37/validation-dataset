#!/bin/bash

#   This script processes a TSV file that contains variant information,
#   extracting coordinates and zygosities for the proband, mother, and father's samples
#   by navigating and extracting information from annotation files (Alamut or VEP) on DNAnexus.
#   Output is written in BED format to be used for downstream analysis with BAMSurgeon.
#
#   Usage: ./find_coords.sh input_file.tsv
#   Input: input_file.csv
#   Outputs: bedfiles/proband.bed; bedfiles/mother.bed; bedfiles/father.bed

set -euo pipefail # avoid silent fails

INPUT_FILE=${1}

# sample details file exists on DNAnexus; lists all samples processed
SAMPLE_DETAILS="Resources:/sample_details/sample_details.tsv"
dx cat "$SAMPLE_DETAILS" > temp/sample_details_cached.tsv

# output directory
BED_DIR="bedfiles"

# ENSURE BEDFILES EXIST (and temp files for processing purposes)
mkdir -p "$BED_DIR" temp
: > temp/temp.txt
: > temp/temp_annot.txt
: > temp/temp_fam.txt
: > log.txt
: > success.txt
: > stderr

for i in proband.bed mother.bed father.bed; do
    # create files if they don't exist
    if [ ! -f "bedfiles/$i" ]; then
    touch bedfiles/$i 
    fi
done


# HELPER FUNCTIONS
extract_zygosity() {
    local gt=$1
    if [[ "$gt" == "0/0" || "$gt" == "0|0" ]]; then
        echo "0"
    elif [[ "$gt" == "0/1" || "$gt" == "1/0" || "$gt" == "0|1" || "$gt" == "1|0" ]]; then
        echo "0.5"
    elif [[ "$gt" == "1/1" || "$gt" == "1|1" ]]; then
        echo "1"
    else
        echo "NA" 
    fi
}

alamut_or_vep() {
    local i="$1"
    alamut=$(dx find data --name "*$i*.annotated.txt" | grep -oP '\(\K[^)]*(?=\))')
    vep=$(dx find data --name "*$i*.vep.vcf.gz" | grep -oP '\(\K[^)]*(?=\))')

    if [[ -n "$alamut" ]]; then echo "alamut"
    elif [[ -n "$vep" ]]; then echo "vep"
    else echo "none"
    fi
}

write_bed() {
    local proband_z=$(extract_zygosity "$PROBAND_GT")
    local mother_z=$(extract_zygosity "$MOTHER_GT")
    local father_z=$(extract_zygosity "$FATHER_GT")

    if [[ "$PROBAND_Z" == "NA" || "$MOTHER_Z" == "NA" || "$FATHER_Z" == "NA" ]]; then
        echo "Skipping row $row_number in $FOLDERNO - zygosity is NA!" >> stderr
        continue
    fi

    echo "Sample info extracted, writing to temp files" >> log.txt

    echo -e "$CHR\t$START\t$END\t$proband_z\t$ALT_ALLELE" >> temp/proband.tmp
    echo -e "$CHR\t$START\t$END\t$mother_z\t$ALT_ALLELE" >> temp/mother.tmp
    echo -e "$CHR\t$START\t$END\t$father_z\t$ALT_ALLELE" >> temp/father.tmp

    sort -u temp/proband.tmp >> "$BED_DIR/proband.bed"
    sort -u temp/mother.tmp >> "$BED_DIR/mother.bed"
    sort -u temp/father.tmp >> "$BED_DIR/father.bed"

    return 0
}


# extract columns from the dataset
tail -n +2 "$INPUT_FILE" | while IFS=$'\t' read -r ROW; do
    FOLDERNO=$(echo "$ROW" | cut -f1)
    CDNACHANGE=$(echo "$ROW" | cut -f3)
    GENENAME=$(echo "$ROW" | cut -f4)
    ALT_ALLELE=$(echo "$ROW" | cut -f8)
    
    printf "\nProcessing folder: %s, %s, %s\n" "$FOLDERNO" "$GENENAME" "$CDNACHANGE" >> log.txt

# get family number from sample details e.g. F08796
FAM_NUM=$(cat temp/sample_details_cached.tsv | grep "$FOLDERNO" | cut -f2)

echo "Family number: $FAM_NUM" >> log.txt


ANNOT_TYPE=$(alamut_or_vep "$FOLDERNO")

# if-else loop choosing between either alamut or vep processing

if [[ "$ANNOT_TYPE" == "alamut" ]]; then

################
#   ALAMUT
################

echo "Processing $ANNOT_TYPE annotation file for $FOLDERNO..." >> log.txt

# Find Alamut file
ANNOTATION_FILE=$(dx find data --name "*$FOLDERNO*.annotated.txt" | grep -oP '\(\K[^)]*(?=\))')
echo "$ANNOTATION_FILE" >> log.txt

# find the relevant columns (matching our gene name and cDNA change) and save to temp.txt
dx cat "$ANNOTATION_FILE" | awk -v gene="$GENENAME" -v cdna="$CDNACHANGE" 'NR==1 || ($0 ~ gene && $0 ~ cdna)' > temp/temp.txt

# get relevant rows: chr.no, gene name, start, end, c.nomen 
tail -n +2 temp/temp.txt | awk -F '\t' '{print $2, $7, $20, $21, $25}' > temp/temp_annot.txt

# FINDING ZYGOSITY INFO 

# get IDs from sample details file
cat temp/sample_details_cached.tsv | grep $FAM_NUM > temp/temp_fam.txt

# Extract sample IDs for proband, mother and father from sample details
while IFS=$'\t' read -r ROW; do
    RELATIONSHIP=$(echo "$ROW" | tr '[:upper:]' '[:lower:]')
    if echo "$RELATIONSHIP" | grep -q "proband"; then
        PROBAND=$(echo "$ROW" | cut -f1)
    elif echo "$RELATIONSHIP" | grep -q "mother"; then
        MOTHER=$(echo "$ROW" | cut -f1)
    elif echo "$RELATIONSHIP" | grep -q "father"; then
        FATHER=$(echo "$ROW" | cut -f1)
    fi
done < temp/temp_fam.txt

HEADER=$(head -n 1 temp/temp.txt)

# convert header line into an array of columns
IFS=$'\t' read -r -a COLUMNS <<< "$HEADER"

# initialise GT column indices to 0 (for duo/singleton handling)
PROBAND_COL=0; MOTHER_COL=0; FATHER_COL=0
GT_PROBAND="GT (${PROBAND})"
GT_MOTHER="GT (${MOTHER})"
GT_FATHER="GT (${FATHER})"

# Find and assign indexes of GT columns
for i in "${!COLUMNS[@]}"; do
    COLNAME="${COLUMNS[$i]}"
    if [ "$COLNAME" = "$GT_PROBAND" ]; then
        PROBAND_COL=$((i+1)) # bash array indices are zero-based but awk uses one-based indexing, so "i+1"
    elif [ "$COLNAME" = "$GT_MOTHER" ]; then
        MOTHER_COL=$((i+1))
    elif [ "$COLNAME" = "$GT_FATHER" ]; then
        FATHER_COL=$((i+1))
    fi
done

echo "GT columns: $PROBAND_COL, $MOTHER_COL, $FATHER_COL (proband, mother, father)" >> log.txt

GENOTYPE_OK=true
row_number=0
# starting from 2nd line (skipping header row), store the entire row as variable "ROW" (no splitting on spaces, tabs or backslashes) and loop through rows.
while IFS= read -r ROW; do
    ((row_number++))
    
    # Extract genotypes from current row. (If missing (column = 0), set to empty string)
    [[ $PROBAND_COL -gt 0 ]] && PROBAND_GT=$(echo "$ROW" | cut -f"$PROBAND_COL") || PROBAND_GT=""
    [[ $MOTHER_COL -gt 0 ]] && MOTHER_GT=$(echo "$ROW" | cut -f"$MOTHER_COL") || MOTHER_GT=""
    [[ $FATHER_COL -gt 0 ]] && FATHER_GT=$(echo "$ROW" | cut -f"$FATHER_COL") || FATHER_GT=""

    echo "GTs: $PROBAND_GT $MOTHER_GT $FATHER_GT" >> log.txt

    # Check for NA/missing genotypes - if so, log and don't write to BED.
    if { [[ $PROBAND_COL -gt 0 && ( -z "$PROBAND_GT" || "$PROBAND_GT" == "NA" ) ]] || \
         [[ $MOTHER_COL -gt 0 && ( -z "$MOTHER_GT" || "$MOTHER_GT" == "NA" ) ]] || \
         [[ $FATHER_COL -gt 0 && ( -z "$FATHER_GT" || "$FATHER_GT" == "NA" ) ]]; }; then
        GENOTYPE_OK=false
        echo "Skipping row $row_number in $FOLDERNO - invalid genotypes!" >> stderr
        continue  # skip to next variant
    fi

    # Extract annotation info for the variant from annot.txt
    ANNOT_ROW=$(sed -n "${row_number}p" temp/temp_annot.txt)
    CHR=$(echo "$ANNOT_ROW" | awk '{print $1}')
    START=$(echo "$ANNOT_ROW" | awk '{print $3}')
    END=$(echo "$ANNOT_ROW" | awk '{print $4}')

    if ! write_bed; then
        continue
    fi

done < <(tail -n +2 temp/temp.txt)


elif [[ "$ANNOT_TYPE" == "vep" ]]; then

################
#   VEP
################

echo "Processing $ANNOT_TYPE annotation file for $FOLDERNO..." >> log.txt

# Find VEP file
ANNOTATION_FILE=$(dx find data --name "*$FOLDERNO*.vep.vcf.gz" | grep -oP '\(\K[^)]*(?=\))')
echo "$ANNOTATION_FILE" >> log.txt

# zcat and wipe metadata lines, then grep for variants that match gene name and DNA change
dx cat "$ANNOTATION_FILE" | zcat | \
    awk -v gene="$GENENAME" -v cdna="$CDNACHANGE" \
    '/^#CHROM/ {print; next} !/^#/ && $0 ~ gene && $0 ~ cdna' > temp/temp.txt

# get relevant rows from temp.txt and pipe into temp_annot.txt
awk -F '\t' -v gene="$GENENAME" -v cdna="$CDNACHANGE" '!/^#/ {print $1, gene, $2, $2, cdna}' temp/temp.txt > temp/temp_annot.txt

# FINDING ZYGOSITY INFO

# get IDs from sample details file
cat temp/sample_details_cached.tsv | grep $FAM_NUM > temp/temp_fam.txt

# Extract sample IDs for proband, mother and father from sample details
while IFS=$'\t' read -r ROW; do
    RELATIONSHIP=$(echo "$ROW" | tr '[:upper:]' '[:lower:]') 

    if echo "$RELATIONSHIP" | grep -q "proband"; then
        PROBAND=$(echo "$ROW" | cut -f1)
    elif echo "$RELATIONSHIP" | grep -q "mother"; then
        MOTHER=$(echo "$ROW" | cut -f1)
    elif echo "$RELATIONSHIP" | grep -q "father"; then
        FATHER=$(echo "$ROW" | cut -f1)
    fi
done < temp/temp_fam.txt

HEADER=$(head -n 1 temp/temp.txt)

# convert header line into an array of columns
IFS=$'\t' read -r -a COLUMNS <<< "$HEADER"

# initialise GT column indices to 0 (for duo/singleton handling)
PROBAND_COL=0;MOTHER_COL=0;FATHER_COL=0

# Find and assign indexes of GT columns
for i in "${!COLUMNS[@]}"; do
    COLNAME="${COLUMNS[$i]}"
    if [ "$COLNAME" = "$PROBAND" ]; then
        PROBAND_COL=$((i+1)) 
    elif [ "$COLNAME" = "$MOTHER" ]; then
        MOTHER_COL=$((i+1))
    elif [ "$COLNAME" = "$FATHER" ]; then
        FATHER_COL=$((i+1))
    fi
done

echo "GT columns: $PROBAND_COL, $MOTHER_COL, $FATHER_COL (proband, mother, father)" >> log.txt

# WRITE INFO TO BED FILES
GENOTYPE_OK=true
row_number=0 
while IFS= read -r ROW; do 
    ((row_number++))

    # Extract genotypes from current row (first sub-field in colon-separated genotype field)
    [[ $PROBAND_COL -gt 0 ]] && PROBAND_GT=$(echo "$ROW" | cut -f"$PROBAND_COL" | cut -d':' -f1) || PROBAND_GT=""
    [[ $MOTHER_COL -gt 0 ]] && MOTHER_GT=$(echo "$ROW" | cut -f"$MOTHER_COL" | cut -d':' -f1) || MOTHER_GT=""
    [[ $FATHER_COL -gt 0 ]] && FATHER_GT=$(echo "$ROW" | cut -f"$FATHER_COL" | cut -d':' -f1) || FATHER_GT=""

    echo "GTs: $PROBAND_GT $MOTHER_GT $FATHER_GT" >> log.txt

    if { [[ $PROBAND_COL -gt 0 && ( -z "$PROBAND_GT" || "$PROBAND_GT" == "NA" ) ]] || \
         [[ $MOTHER_COL -gt 0 && ( -z "$MOTHER_GT" || "$MOTHER_GT" == "NA" ) ]] || \
         [[ $FATHER_COL -gt 0 && ( -z "$FATHER_GT" || "$FATHER_GT" == "NA" ) ]]; }; then
        GENOTYPE_OK=false
        echo "Skipping row $row_number in $FOLDERNO - invalid genotypes!" >> stderr
        continue  # skip to next variant
    fi

    # Extract annotation info for the variant from annot.txt
    ANNOT_ROW=$(sed -n "${row_number}p" temp/temp_annot.txt)
    CHR=$(echo "$ANNOT_ROW" | awk '{print $1}')
    START=$(echo "$ANNOT_ROW" | awk '{print $3}')
    END=$(echo "$ANNOT_ROW" | awk '{print $4}')

    if ! write_bed; then
        continue
    fi

done < <(tail -n +2 temp/temp.txt)

else
    echo "ERROR: Could not detect annotation type for $FOLDERNO" >> stderr
    continue
fi

# clean up logs
MIN_EXPECTED_LINES=8 # minimum lines that would be printed per successful run (accounting for duplicates)
ACTUAL_LINES=$(wc -l < log.txt)

# if script fails partway through processing a sample, send this to stderr
if [ "$ACTUAL_LINES" -ge "$MIN_EXPECTED_LINES" ]; then
    echo "$FOLDERNO" >> success.txt
else
    cat log.txt >> stderr
fi

> log.txt # wipe log file for the next sample

done