# Validation Dataset

> This repository is under active development

The goal of this project is to create a validation cohort using [BamSurgeon](https://github.com/adamewing/bamsurgeon), a tool that allows variants to be inserted into preexisting BAM files. By creating "FrankenBAM" samples with spliced-in variants, we can validate hundreds of variants within a single sample, significantly reducing both the time and computational resources required for validation purposes. 

---
## Prerequisites
- Bash shell
- Docker image of BamSurgeon (`addsnv.py`)
---

## 1. Obtain variant list

A list of pathogenic/likely pathogenic variants previously identified in the laboratory and their metadata was generated using Crystal Reports. Resulting file is a `csv`.

## 2. Convert dataset into format compatible with BamSurgeon

In order to use these variants as our BamSurgeon input, they must be formatted as a BED file with the following columns: 
`chr    start   end VAF alt`

Preprocessing scripts (`/preprocess`):

R script `modifydata.R`:
- Cleans Crystal Report `csv`, removing duplicated/incomplete entries
- Extracts SNVs (removing indels which are not accepted by BamSurgeon's `addsnv.py`)
- Outputs new `tsv` files with sample lookup details for DNAnexus and variant metadata
- Produces separate male and female `tsv` files to allow for creation of two family trios (for accurate handling of X-linked variants)

Bash script `annotation_to_bed.sh`:
- Navigates through cloud server to locate annotation files for old samples
- Script pulls relevant info (chromosome, start coord, end coord, zygosity)
    - Annotation files for samples will have either been created by VEP or Alamut, depending on when the sample was analysed. VEP and Alamut annotation files have unique formats, so are handled differently. 

As exomes/genomes are run as family trios, we also require a BED file representing the mother and father. 
- Logic for duos: the other parent is effectively assumed to match the reference genome, since the variant will not be added to the corresponding BED file. 

## 3. Running BamSurgeon

Using BamSurgeon's `addsnv.py` script to splice in SNVs. 

Command for bamsurgeon:
```bash
python3 /bamsurgeon/bin/addsnv.py  -v input.bed -f input.bam --aligner mem --picardjar /picard.jar -p 8 -o output.bam -r ref.fasta
```

#### WORK IN PROGRESS:
Currently, our aim is to implement and test a workflow for running `addsnv.py` across four families of trios (one male proband and one female proband, for both WGS and WES data).
- The bash script implementation `process_bams.sh` was a first attempt that loops through samples, but may be too slow for large datasets. 
- A WDL workflow (`/WDL`) is being tested as a more scalable implementation (in development).

---

Project being undertaken as a component of my STP training. 
