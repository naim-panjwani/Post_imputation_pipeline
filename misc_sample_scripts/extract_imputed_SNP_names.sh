#!/bin/bash
touch 17-1-imputed_SNP_list.txt
for i in {1..22}; do
  tabix 16-JME_Round1_and_2_chr${i}_beagle5_imputed.vcf.gz ${i}: |cut -f3 >>17-1-imputed_SNP_list.txt
done
i="X"
tabix 16-JME_Round1_and_2_chr${i}_beagle5_imputed.vcf.gz ${i}: |cut -f3 >>17-1-imputed_SNP_list.txt
i="Y"
tabix 16-JME_Round1_and_2_chr${i}_beagle5_imputed.vcf.gz ${i}: |cut -f3 >>17-1-imputed_SNP_list.txt
i="MT"
tabix 16-JME_Round1_and_2_chr${i}_beagle5_imputed.vcf.gz ${i}: |cut -f3 >>17-1-imputed_SNP_list.txt
