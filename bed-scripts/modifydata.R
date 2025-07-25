#!/usr/bin/env Rscript

# modifydata.R
# This script takes a CSV of pathogenic or LP variants from a Crystal Report
# and cleans it. The output TSV can be used as an input for annotation_to_bed.sh. 

#-------------------------------------------
# LOAD PACKAGES AND DATA
#-------------------------------------------

library(tidyverse)

# import data
data <- R14_WGS_variants_for_validation.csv

# check data
glimpse(data)
head(data)

#-------------------------------------------
# CLEAN DATASET
#-------------------------------------------

### reformat/remove duplicates
data <- data %>%
  distinct() %>%
  mutate(
    GENE = str_to_upper(GENE) %>% # convert gene names to uppercase
    str_replace_all(" ", "_"), # replace spaces with underscores
  CDNACHANGE = as.character(CDNACHANGE), 
  ZYGOSITY = as.character(ZYGOSITY)
  )


#--------------------------------------------
# EXTRACT ALT ALLELE
#--------------------------------------------

alt_allele <- function(cdna_change) {
  match <- stringr::str_match(cdna_change, "^c\\.[0-9_]+[ACGT]>[ACGT]$")
  if (is.na(match[1])) return(NA_character_)
  
  parts <- stringr::str_match(cdna_change, "^c\\.[0-9_]+([ACGT])>([ACGT])$")
  return(parts[,3]) 
}

#--------------------------------------------
# CREATE ZYGOSITY FIELD 
#--------------------------------------------

data2 <- data %>%
  mutate(
    AF = case_when(
      str_to_lower(ZYGOSITY) == "heterozygous" ~ 0.5,
      str_to_lower(ZYGOSITY) == "homozygous" ~ 1.0,
      TRUE ~ NA_real_
    ),
    ALT_ALLELE = alt_allele(CDNACHANGE)
  )

#--------------------------------------------
# EXTRACT RELEVANT FIELDS FOR GENE COORDINATE LOOKUPS
#--------------------------------------------

data3 <- data2 %>%
  select(
    FOLDERNO,
    ORDNO,
    CDNACHANGE,
    GENE,
    GENDER,
    ZYGOSITY,
    RELATIONSHIP,
    ALT_ALLELE
  ) %>%
  drop_na(ALT_ALLELE)


#--------------------------------------------
# SPLIT INTO FAMILIES
#--------------------------------------------
 
# Male probands
male_proband <- data3 %>%
  filter(
    str_to_lower(RELATIONSHIP) == "proband",
    str_to_lower(GENDER)       == "male"
  )

# Female probands
female_proband <- data3 %>%
  filter(
    str_to_lower(RELATIONSHIP) == "proband",
    str_to_lower(GENDER)       == "female"
  )

# check tables:
list(
  male_proband   = glimpse(male_proband),
  female_proband = glimpse(female_proband),
)

#--------------------------------------------
# RUN BASH SCRIPT TO GET COORDS
#--------------------------------------------
# The bash script will create three bed files: proband.bed, mother.bed, father.bed.

# Create output dirs
dir.create("output")
dir.create("output/male")
dir.create("output/female")

# save tables as TSVs
write_tsv(male_proband, "output/male/input.tsv")
write_tsv(female_proband, "output/female/input.tsv")

# Run script for male proband
message("Running annotation_to_bed.sh for male proband...")
system("bash annotation_to_bed.sh output/male/input.tsv output/male")

# Run script for female proband
message("Running annotation_to_bed.sh for female proband...")
system("bash annotation_to_bed.sh output/female/input.tsv output/female")
