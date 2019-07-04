#!/usr/bin/R
# Script to check whether the SNPs that we want to add back are SNPs that are not already present in the imputed files
# These "undesired SNPs" will be listed out and won't be extracted from the genotypes
# Syntax: Rscript 01-cleanup_SNPs_to_add_back.R <SNP_list_to_add_back_filename> <imputed_SNP_list_filename> <output_filename>

library(argparser)
library(data.table)

p <- arg_parser("Cleanup tool to check whether the SNPs that we want to add back are SNPs that are not already present in the imputed files") 

p <- add_argument(p, "SNP_list_to_add_back_filename", help="Text file with the list of SNPs we would like to add back")
p <- add_argument(p, "imputed_SNP_list_filename", help="A one-column text file of all SNP ID's in the imputed files")
p <- add_argument(p, "output_filename", help="Desired output filename")

argv <- parse_args(p)
SNP_list_to_add_back_filename <- argv$SNP_list_to_add_back_filename
imputed_SNP_list_filename <- argv$imputed_SNP_list_filename
output_filename <- argv$output_filename

#SNP_list_to_add_back_filename <- "17-Omni25_SNPs_to_add_back.txt"
#imputed_SNP_list_filename <-"17-1-imputed_SNP_list.txt"
#output_filename <- "17-2-Omni25_SNPs_to_add_back.txt"

print("Reading SNP list files")
SNP_list_to_add_back <- fread(SNP_list_to_add_back_filename, header=F, stringsAsFactors=F)
imputed_SNP_list <- fread(imputed_SNP_list_filename, header=F, stringsAsFactors=F)

print("Listing final set of genotyped SNPs in the genotyped SNP list that are not imputed")
not_imputed <- SNP_list_to_add_back[-which(SNP_list_to_add_back$V1 %in% imputed_SNP_list$V1),]
num_not_imputed <- nrow(not_imputed)
print(paste0("Found ", num_not_imputed, " genotyped SNPs that should be added back"))

print(paste0("Saving list of genotyped SNPs to be added back to the imputed files in ", output_filename))
fwrite(not_imputed, output_filename, quote=F, row.names=F, col.names=F)


