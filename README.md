# Post imputation pipeline
This pipeline aims to add genotyped SNPs that get thrown out due to their absence from the reference back into the imputed VCFs post-imputation (specifically, post-imputation with BEAGLE v4 or later)

## Step 0
You will need the following files, programs, libraries:
- File containing the list of SNPs to add back (can be raw)
- File containing the list of SNPs that are imputed (check out the misc\_sample\_scripts/extract\_imputed\_SNP\_names.sh sample script)
- Will Rayner's strand alignment file for your specific platform (see https://www.well.ox.ac.uk/~wrayner/strand/)
- The 1000 Genomes reference fasta file
- R packages: argparse, data.table
- bgzip, tabix, and vcftools

## Step 1
- Check list of SNPs to add back to ensure that they are not already imputed
```
Rscript 01-cleanup\_SNPs\_to\_add\_back.R <SNP\_list\_to\_add\_back\_filename> <imputed\_SNP\_list\_filename> <output\_filename>
Rscript Post_imputation_pipeline/01-cleanup_SNPs_to_add_back.R 17-Omni25_SNPs_to_add_back.txt 17-1-imputed_SNP_list.txt 17-2-Omni25_SNPs_to_add_back.txt
```
For the SNP list to add back, you could simply list out all SNPs in the bim file prior to the conform-gt step

## Step 2
- Extract the SNPs using PLINK
plink --bfile 08-1-underscore\_to\_dash\_id\_update --extract 17-2-Omni25\_SNPs\_to\_add\_back.txt --make-bed --out 17-SNPs\_to\_add\_back


