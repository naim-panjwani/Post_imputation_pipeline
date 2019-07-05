# Post imputation pipeline
This pipeline aims to add genotyped SNPs that get thrown out due to their absence from the reference back into the imputed VCFs post-imputation (specifically, post-imputation with BEAGLE v4 or later). BEAGLE version 5 was used here.

## Step 0
You will need the following files, programs, libraries:
- File containing the list of SNPs to add back (can be raw)
- File containing the list of SNPs that are imputed (check out the misc\_sample\_scripts/extract\_imputed\_SNP\_names.sh sample script)
- Will Rayner's strand alignment file for your specific platform (see https://www.well.ox.ac.uk/~wrayner/strand/)
- The 1000 Genomes reference fasta file
- R packages: argparse, data.table
- bgzip, tabix, vcftools, and bcftools
- PLINK (v1.90b6.9 64-bit (4 Mar 2019) used here)

## Step 1
- Check list of SNPs to add back to ensure that they are not already imputed
```
Rscript 01-cleanup_SNPs_to_add_back.R <SNP_list_to_add_back_filename> <imputed_SNP_list_filename> <output_filename>
```
Example:
```
Rscript Post_imputation_pipeline/01-cleanup_SNPs_to_add_back.R 17-Omni25_SNPs_to_add_back.txt 17-1-imputed_SNP_list.txt 17-2-Omni25_SNPs_to_add_back.txt
```
For the SNP list to add back, you could simply list out all SNPs in the bim file prior to the conform-gt step.  
You may invoke help on the script via the -h argument:
```
Rscript 01-cleanup_SNPs_to_add_back.R -h
usage: 01-cleanup_SNPs_to_add_back.R [--] [--help] [--opts OPTS] SNP_list_to_add_back_filename imputed_SNP_list_filename output_filename

Cleanup tool to check whether the SNPs that we want to add back are SNPs that are not already present in the imputed files

positional arguments:
  SNP_list_to_add_back_filename			Text file with the list of SNPs we would like to add back
  imputed_SNP_list_filename			A one-column text file of all SNP ID's in the imputed files
  output_filename			Desired output filename

flags:
  -h, --help			show this help message and exit

optional arguments:
  -x, --opts OPTS			RDS file containing argument values
```


## Step 2
- Extract the SNPs using PLINK
```
plink --bfile <bplink_filename> --extract <step1_snp_list_output> --make-bed --out <output_filename>
```
Example:
```
plink --bfile 08-1-underscore_to_dash_id_update --extract 17-2-Omni25_SNPs_to_add_back.txt --make-bed --out 17-SNPs_to_add_back
```


## Step 3
- Align the SNPs to the plus strand
Since the Illumina SNP names were converted to rsid's, we need to convert them back for strand alignment using Will Rayner's update\_build.sh script
```
bash 03-align_strand.sh <bPLINK_filename> <rsid_illumina_mapping_filename> <strand_filename> <output_filename>
```
Example:
```
bash 03-align_strand.sh 17-SNPs_to_add_back 18-SNP_renaming_1-3_and_1-4_merged.txt strand_files_1-4/InfiniumOmni2-5-8v1-4_A1-b37.strand 20-Strand_aligned_1-4
```

## Step 4
- Convert to VCF; fix the X,Y,MT and pseudoautosomal X chrom field; and align the reference allele to 1000 Genomes
```
bash 04-convert_and_fix_vcf.sh <bPLINK_filename> <output_filename>
```
Example:
```
bash 04-convert_and_fix_vcf.sh 20-Strand_aligned_1-4 23-SNPs_to_add_back_normalized.vcf.gz
```

## Step 5
- Format the VCF to match the format of BEAGLE (version 5) output VCF
```
python3 05-VCF_Reformatting.py <VCF_filename> <output_filename>
```
Example:
```
python3 Post_imputation_pipeline/05-VCF_Reformating.py 23-SNPs_to_add_back_normalized.vcf.gz 24-SNPs_to_add_back_reformatted.vcf
```
Script help details:
```
python3 05-VCF_Reformating.py -h
usage: 05-VCF_Reformating.py [-h] vcf_filename output_filename

Tool to format a raw VCF to match the format of BEAGLE (v5) output VCF

positional arguments:
  vcf_filename     The VCF filename that needs reformatting
  output_filename  Output filename

optional arguments:
  -h, --help       show this help message and exit
```

## Step 6
- Fix, sort and compress
```
sed 's/^\t//g' <step5_outputfilename> |sort -k1,1V -k2,2n |bgzip -c ><step5_outputfilename>.gz
tabix -p vcf <step5_outputfilename>.gz
if [ -e <step5_outputfilename>.gz ]; then rm <step5_outputfilename>; fi
```
Example:
```
sed 's/^\t//g' 24-SNPs_to_add_back_reformatted.vcf |sort -k1,1V -k2,2n |bgzip -c >24-SNPs_to_add_back_reformatted.vcf.gz
tabix -p vcf 24-SNPs_to_add_back_reformatted.vcf.gz
if [ -e 24-SNPs_to_add_back_reformatted.vcf.gz.tbi ]; then rm 24-SNPs_to_add_back_reformatted.vcf; fi
```

## Step 7
- Merge and sort with the imputed files. If the imputed files are split by chromosome, then an example of how to combine these follows:
```
for i in {1..22}; do
  (tabix -h 16-JME_Round1_and_2_chr${i}_beagle5_imputed.vcf.gz ${i}: ; tabix 24-SNPs_to_add_back_reformatted.vcf.gz ${i}: ) |vcf-sort |bgzip -c >25-JME_Round1_and_2_chr${i}_imputed_all_snps_in.vcf.gz
  tabix -p vcf 25-JME_Round1_and_2_chr${i}_imputed_all_snps_in.vcf.gz
done
i="X"
(tabix -h 16-JME_Round1_and_2_chr${i}_beagle5_imputed.vcf.gz ${i}: ; tabix 24-SNPs_to_add_back_reformatted.vcf.gz ${i}: ) |vcf-sort |bgzip -c >25-JME_Round1_and_2_chr${i}_imputed_all_snps_in.vcf.gz
25-JME_Round1_and_2_chr${i}_imputed_all_snps_in.vcf.gz
i="Y"
tabix -h 24-SNPs_to_add_back_reformatted.vcf.gz ${i}: |vcf-sort |bgzip -c >25-JME_Round1_and_2_chr${i}_imputed_all_snps_in.vcf.gz
25-JME_Round1_and_2_chr${i}_imputed_all_snps_in.vcf.gz
i="MT"
tabix -h 24-SNPs_to_add_back_reformatted.vcf.gz ${i}: |vcf-sort |bgzip -c >25-JME_Round1_and_2_chr${i}_imputed_all_snps_in.vcf.gz
25-JME_Round1_and_2_chr${i}_imputed_all_snps_in.vcf.gz
```
The 16-JME\_Round1\_and\_2\_chr${i}\_beagle5\_imputed.vcf.gz are the imputed files, and 24-SNPs\_to\_add\_back\_reformatted.vcf.gz is the file generated from step 6

