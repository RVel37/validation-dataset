#!/bin/bash

INPUT_FILE=${1}
SAMPLE_DETAILS="Resources:/sample_details/sample_details.tsv"
BED_DIR="bedfiles"

# ENSURE BEDFILES EXIST
mkdir -p bedfiles; mkdir -p temp
touch temp/temp.txt; touch temp/temp_annot.txt; touch temp/temp_fam.txt

for i in proband.bed mother.bed father.bed; do
    if [ ! -f "bedfiles/$i" ]; then
    touch bedfiles/$i # create files if they don't exist
    exit 1
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


# extract useful columns from input csv
tail -n +2 "$INPUT_FILE" | while IFS=$'\t' read -r ROW; do
    FOLDERNO=$(echo "$ROW" | cut -f1)
    CDNACHANGE=$(echo "$ROW" | cut -f3)
    GENENAME=$(echo "$ROW" | cut -f4)

echo "Processing folder: $FOLDERNO, $GENENAME, $CDNACHANGE"

# get family number from sample details e.g. F08796
FAM_NUM=$(dx cat $SAMPLE_DETAILS | grep $FOLDERNO | cut -f2)

echo "Family number: $FAM_NUM"

# CHECK WHETHER ALAMUT OR VEP ANNOTATION FILE



# get annotation file (in brackets)
ANNOTATION_FILE=$(dx find data --name "*$FOLDERNO*.annotated.txt" | grep -oP '\(\K[^)]*(?=\))')

# find relevant columns and save
dx cat "$ANNOTATION_FILE" | awk -v gene="$GENENAME" -v cdna="$CDNACHANGE" 'NR==1 || ($0 ~ gene && $0 ~ cdna)' > temp/temp.txt

# get relevant rows: chr.no, gene name, start, end, c.nomen 
awk -F '\t' '{print $2, $7, $20, $21, $25}' temp/temp.txt > temp/temp_annot.txt 

# get zygosities of parents
# from temp.txt pull out all columns GT(*) 

# get mother and father's IDs
dx cat Resources:/sample_details/sample_details.tsv | grep $FAM_NUM > temp/temp_fam.txt

PROBAND=$(awk '$3=="Proband" {print $1}' "temp/temp_fam.txt")
MOTHER=$(awk '$NF=="mother" {print $1}' "temp/temp_fam.txt")
FATHER=$(awk '$NF=="father" {print $1}' "temp/temp_fam.txt")

# match each GT column with PROBAND, MOTHER and FATHER 

# GT column names:
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

echo "Sample IDs:"
echo "Proband: $PROBAND"
echo "Mother:  $MOTHER"
echo "Father:  $FATHER"


# header line
HEADER=$(head -n 1 temp/temp.txt)

# create an array from the header columns + loop through to find GT col numbers
IFS=$'\t' read -r -a COLUMNS <<< "$HEADER"

GT_PROBAND="GT (${PROBAND})"
GT_MOTHER="GT (${MOTHER})"
GT_FATHER="GT (${FATHER})"

for i in "${!COLUMNS[@]}"; do
    COLNAME="${COLUMNS[$i]}"
    if [ "$COLNAME" = "$GT_PROBAND" ]; then
        PROBAND_COL=$((i+1))
    elif [ "$COLNAME" = "$GT_MOTHER" ]; then
        MOTHER_COL=$((i+1))
    elif [ "$COLNAME" = "$GT_FATHER" ]; then
        FATHER_COL=$((i+1))
    fi
done

echo "GT columns: $PROBAND_COL, $MOTHER_COL, $FATHER_COL (proband, mother, father)"

awk -F'\t' -v p="$PROBAND_COL" -v m="$MOTHER_COL" -v f="$FATHER_COL" '{print $p, $m, $f}' temp/temp.txt


row_number=0
tail -n +2 temp/temp.txt | while IFS= read -r ROW; do
        ((row_number++))

        # Extract GTs from the current row
        PROBAND_GT=$(echo "$ROW" | cut -f"$PROBAND_COL")
        MOTHER_GT=$(echo "$ROW" | cut -f"$MOTHER_COL")
        FATHER_GT=$(echo "$ROW" | cut -f"$FATHER_COL")

        # Debugging
        echo "GTs: $PROBAND_GT $MOTHER_GT $FATHER_GT"

        # Get annotation for this row
        ANNOT_ROW=$(sed -n "${row_number}p" temp/temp_annot.txt)
        CHR=$(echo "$ANNOT_ROW" | awk '{print $1}')
        START=$(echo "$ANNOT_ROW" | awk '{print $3}')
        END=$(echo "$ANNOT_ROW" | awk '{print $4}')

        # Convert to zygosity
        PROBAND_Z=$(extract_zygosity "$PROBAND_GT")
        MOTHER_Z=$(extract_zygosity "$MOTHER_GT")
        FATHER_Z=$(extract_zygosity "$FATHER_GT")

        # Write to BED files 

        ### ERROR: THIS PRINTS THE HEADER TOO. Do tail -n -2 ? ###
        echo -e "$CHR\t$START\t$END\t$PROBAND_Z" >> "$BED_DIR/proband.bed"
        echo -e "$CHR\t$START\t$END\t$MOTHER_Z" >> "$BED_DIR/mother.bed"
        echo -e "$CHR\t$START\t$END\t$FATHER_Z" >> "$BED_DIR/father.bed"
    
    ### NEED FUNCTION TO PIPE FAILURES (e.g. quads/duos??) ###
    ### for duos, skipping the parent means they are assumed to match the reference genome. 

    done

done
