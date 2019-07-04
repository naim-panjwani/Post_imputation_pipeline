# Post imputation pipeline
This pipeline aims to add genotyped SNPs that get thrown out due to their absence from the reference back into the imputed VCFs post-imputation (specifically, post-imputation with BEAGLE v4 or later)

## Step 1
- Check list of SNPs to add back to ensure that they are not already imputed
```
Rscript 01-cleanup_SNPs_to_add_back.R <SNP_list_to_add_back_filename> <imputed_SNP_list_filename> <output_filename>
```

