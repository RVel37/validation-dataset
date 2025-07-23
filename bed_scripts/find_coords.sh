#!/bin/bash

#   This script processes a CSV file that contains variant information,
#   extracting coordinates and zygosities for the proband, mother, and father's samples
#   by navigating and extracting information from annotation files (Alamut or VEP) on DNAnexus.
#   Output is written in BED format to be used for downstream analysis with BAMSurgeon.
#
#   Usage: ./find_coords.sh input_file.csv
#   Input: input_file.csv
#   Outputs: bedfiles/proband.bed; bedfiles/mother.bed; bedfiles/father.bed


INPUT_FILE=${1} 
# sample details file exists on DNAnexus; lists all samples processed
SAMPLE_DETAILS="Resources:/sample_details/sample_details.tsv"
# output directory
BED_DIR="bedfiles"

# ENSURE BEDFILES EXIST (and temp files for processing purposes)
mkdir -p bedfiles; mkdir -p temp
touch temp/temp.txt; touch temp/temp_annot.txt; touch temp/temp_fam.txt

for i in proband.bed mother.bed father.bed; do
    if [ ! -f "bedfiles/$i" ]; then
    touch bedfiles/$i # create files if they don't exist
    fi
done


# Helper function
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


# extract columns from the dataset
tail -n +2 "$INPUT_FILE" | while IFS=$'\t' read -r ROW; do
    FOLDERNO=$(echo "$ROW" | cut -f1)
    CDNACHANGE=$(echo "$ROW" | cut -f3)
    GENENAME=$(echo "$ROW" | cut -f4)

echo "Processing folder: $FOLDERNO, $GENENAME, $CDNACHANGE" >> log.txt

# get family number from sample details e.g. F08796
FAM_NUM=$(dx cat "$SAMPLE_DETAILS" | grep "$FOLDERNO" | cut -f2)

echo "Family number: $FAM_NUM" >> log.txt

# CHECK WHETHER ALAMUT OR VEP ANNOTATION FILE
alamut_or_vep() {
    local i="$1"
    alamut=$(dx find data --name "*$i*.annotated.txt" | grep -oP '\(\K[^)]*(?=\))')
    vep=$(dx find data --name "*$i*.vep.vcf.gz" | grep -oP '\(\K[^)]*(?=\))')

    if [[ -n "$alamut" ]]; then
        echo "alamut"
    elif [[ -n "$vep" ]]; then
        echo "vep"
    else
        echo "none"
    fi
}

ANNOT_TYPE=$(alamut_or_vep "$FOLDERNO")
echo "annotation type = $ANNOT_TYPE"

################
#   VEP
################

echo "Using VEP annotation file for $FOLDERNO" >> log.txt

# Find VEP file (zipped vcf)
VCF=$(zcat "*$FOLDERNO*.vep.vcf.gz")
ANNOTATION_FILE=$(dx find data --name "*$VCF" | grep -oP '\(\K[^)]*(?=\))')

# zcat and wipe metadata lines, then grep for variants that match gene name and DNA change
dx cat "$ANNOTATION_FILE" | zcat | grep -v "^##" | grep "$GENENAME" | grep "$CDNACHNAGE" > temp/temp.txt

# get relevant rows from temp.txt: CROM and POS
awk -F '\t' -v gene="$GENENAME" -v cdna="$CDNACHANGE" '{print $1, gene, $2, $2, $2, cdna}'

# FINDING ZYGOSITY INFO
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

# debugging
echo "Sample IDs:"
echo "Proband: $PROBAND"
echo "Mother:  $MOTHER"
echo "Father:  $FATHER"

# header line
HEADER=$(head -n 1 temp/temp.txt)

# convert header line into an array of columns
IFS=$'\t' read -r -a COLUMNS <<< "$HEADER"

# initialise GT column indices to 0 (if not found, i.e. for duos/singletons)
PROBAND_COL=0
MOTHER_COL=0
FATHER_COL=0

# GT header column names
GT_PROBAND="${PROBAND}"
GT_MOTHER="${MOTHER}"
GT_FATHER="${FATHER}"

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

#debugging
echo "GT columns: $PROBAND_COL, $MOTHER_COL, $FATHER_COL (proband, mother, father)"

# pull the GT field from each sample column - if they are 0 then print nothing
awk -F'\t' -v p="$PROBAND_COL" -v m="$MOTHER_COL" -v f="$FATHER_COL" '
{
    if (p>0) {split($p, p, ":"); proband_gt=p[1]}
    else {proband_gt = ""}

    if (m>0) {split($m, m, ":"); mother_gt=m[1]}
    else {mother_gt = ""}

    if (f>0) {split($f, f, ":"); father_gt=f[1]}
    else {father_gt = ""}

    print proband_gt, mother_gt, father_gt
}'  temp/temp.txt


# WRITE INFO TO BED FILES
row_number=0
tail -n +2 temp/temp.txt | while IFS= read -r ROW; do   # starting from 2nd line (skipping header row), store the entire row as variable "ROW" (not splitting on spaces, tabs or backslashes) and loop through





else
    echo "ERROR: Could not find annotation file for $FOLDERNO" >&2 # send to stderr
    continue
fi