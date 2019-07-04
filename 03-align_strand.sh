#!/bin/bash
# Since the Illumina SNP names were converted to rsid's, we need to convert them back for strand alignment using Will Rayner's update_build.sh script
# Syntax: bash 03-align_strand.sh <bPLINK_filename> <rsid_illumina_mapping_filename> <strand_filename> <output_filename>
# Example: bash 03-align_strand.sh 17-SNPs_to_add_back 18-SNP_renaming_1-3_and_1-4_merged.txt strand_files_1-4/InfiniumOmni2-5-8v1-4_A1-b37.strand 20-Strand_aligned_1-4

# 

bplink_filename=$1
rsid_illumina_mapping_filename=$2
strand_filename=$3
output_filename=$4

temp_id="201907041804"

echo "Switching to Illumina SNP ID's"
plink --bfile $bplink_filename --update-name $rsid_illumina_mapping_filename --make-bed --out ${temp_id}-temp1

echo "Aligning to the plus strand"
bash update_build.sh ${temp_id}-temp1 $strand_filename ${temp_id}-temp2

echo "Switching back to rsid's"
awk '{print $2"\t"$1}' $rsid_illumina_mapping_filename >${temp_id}-temp3
plink --bfile ${temp_id}-temp2 --update-name ${temp_id}-temp3 --make-bed --out $output_filename

echo "Removing temporary files"
rm ${temp_id}-temp*

