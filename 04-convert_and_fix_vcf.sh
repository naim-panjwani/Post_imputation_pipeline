#!/bin/bash
# Converts to VCF; fixes the X,Y,MT and pseudoautosomal X chrom field; and align the reference allele to 1000 Genomes
# Input: a binary PLINK file
# Output: a left-aligned and 1KG-normalized VCF
# Dependencies: ensure the 1000 Genomes reference fasta files are present:
#   human_g1k_v37.fasta
#   human_g1k_v37.fasta.fai
#   human_g1k_v37.dict
# Syntax: bash 04-convert_and_fix_vcf.sh <bPLINK_filename> <output_filename>
# Example: bash 04-convert_and_fix_vcf.sh 20-Strand_aligned_1-4 23-SNPs_to_add_back_normalized.vcf.gz

bplink_filename=$1
output_filename=$2

temp_id="201907041825"

echo "Extracting major alleles"
cut -f5,2 ${bplink_filename}.bim >${temp_id}-1
plink --bfile $bplink_filename --recode vcf --a1-allele ${temp_id}-1 --out ${temp_id}-2

echo "Fixing the chromosome field for X, Y, pseudo-X and MT"
sed 's/^23/X/g' ${temp_id}-2.vcf |sed 's/^24/Y/g' |sed 's/^25/X/g' |sed 's/^26/MT/g' |vcf-sort |bgzip -c >${temp_id}-3.vcf.gz
tabix -p vcf ${temp_id}-3.vcf.gz

echo "Aligning the reference allele to 1KG"
bcftools norm --check-ref ws -f human_g1k_v37.fasta ${temp_id}-3.vcf.gz -o $output_filename -O z
tabix -p vcf $output_filename

echo "Removing temporary files"
rm ${temp_id}*
