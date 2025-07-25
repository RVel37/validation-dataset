#!/bin/bash

#   This script processes a TSV file that contains variant information,
#   extracting coordinates and zygosities for the proband, mother, and father's samples
#   by navigating and extracting information from annotation files (Alamut or VEP) on DNAnexus.
#   Output is written in BED format to be used for downstream analysis with BAMSurgeon.
#
#   Usage: ./find_coords.sh input_file.tsv
#   Input: input_file.csv
#   Outputs: bedfiles/proband.bed; bedfiles/mother.bed; bedfiles/father.bed


INPUT_FILE=${1} 
# sample details file exists on DNAnexus; lists all samples processed
SAMPLE_DETAILS="Resources:/sample_details/sample_details.tsv"
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


# extract columns from the dataset
tail -n +2 "$INPUT_FILE" | while IFS=$'\t' read -r ROW; do
    FOLDERNO=$(echo "$ROW" | cut -f1)
    CDNACHANGE=$(echo "$ROW" | cut -f3)
    GENENAME=$(echo "$ROW" | cut -f4)
    
    printf "\nProcessing folder: %s, %s, %s\n" "$FOLDERNO" "$GENENAME" "$CDNACHANGE" >> log.txt

# get family number from sample details e.g. F08796
FAM_NUM=$(dx cat "$SAMPLE_DETAILS" | grep "$FOLDERNO" | cut -f2)

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
dx cat "$SAMPLE_DETAILS" | grep $FAM_NUM > temp/temp_fam.txt

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
tail -n +2 temp/temp.txt | while IFS= read -r ROW; do
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

    # Convert to zygosity - run predefined extract_zygosity function
    PROBAND_Z=$(extract_zygosity "$PROBAND_GT")
    MOTHER_Z=$(extract_zygosity "$MOTHER_GT")
    FATHER_Z=$(extract_zygosity "$FATHER_GT")

    echo "Sample info extracted, writing to BEDfile" >> log.txt

    # Write to BED files
    echo -e "$CHR\t$START\t$END\t$PROBAND_Z" >> "$BED_DIR/proband.bed"
    echo -e "$CHR\t$START\t$END\t$MOTHER_Z" >> "$BED_DIR/mother.bed"
    echo -e "$CHR\t$START\t$END\t$FATHER_Z" >> "$BED_DIR/father.bed"
done


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
dx cat "$SAMPLE_DETAILS" | grep $FAM_NUM > temp/temp_fam.txt

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

# extract genotypes from GT field in each sample column and assign
# this is done by splitting each sample column into a temporary array divided by colons
# read PROBAND_GT MOTHER_GT FATHER_GT < <(
#     awk -F'\t' -v p="$PROBAND_COL" -v m="$MOTHER_COL" -v f="$FATHER_COL" '
#     NR==2 {
#         if (p>0) {split($p, p_split, ":"); proband_gt = p_split[1]}
#         else {proband_gt = ""}

#         if (m>0) {split($m, m_split, ":"); mother_gt = m_split[1]}
#         else {mother_gt = ""}

#         if (f>0) {split($f, f_split, ":"); father_gt = f_split[1]}
#         else {father_gt = ""}

#         print proband_gt, mother_gt, father_gt
#     }' temp/temp.txt
# )

# WRITE INFO TO BED FILES
GENOTYPE_OK=true
row_number=0 
tail -n +2 temp/temp.txt | while IFS= read -r ROW; do 
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

    # Convert to zygosity - run predefined extract_zygosity function
    PROBAND_Z=$(extract_zygosity "$PROBAND_GT")
    MOTHER_Z=$(extract_zygosity "$MOTHER_GT")
    FATHER_Z=$(extract_zygosity "$FATHER_GT")

    echo "Sample info extracted, writing to BEDfile" >> log.txt

    # Write to BED files
    echo -e "$CHR\t$START\t$END\t$PROBAND_Z" >> "$BED_DIR/proband.bed"
    echo -e "$CHR\t$START\t$END\t$MOTHER_Z" >> "$BED_DIR/mother.bed"
    echo -e "$CHR\t$START\t$END\t$FATHER_Z" >> "$BED_DIR/father.bed"
done

else
    echo "ERROR: Could not detect annotation type for $FOLDERNO" >> stderr
    continue
fi

# clean up logs
MIN_EXPECTED_LINES=8 # minimum lines that would be printed per successful run (accounting for duplicates in annotation file)
ACTUAL_LINES=$(wc -l < log.txt)

# if script fails partway through processing a sample, send this to stderr
if [ "$ACTUAL_LINES" -ge "$MIN_EXPECTED_LINES" ]; then
    echo "$FOLDERNO" >> success.txt
else
    cat log.txt >> stderr
fi

> log.txt # wipe log file for the next sample

done