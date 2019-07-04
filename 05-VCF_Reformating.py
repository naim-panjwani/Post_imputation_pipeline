#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 30 11:37:28 2019

@author: naim
"""
import pandas as pd
import numpy as np
import gzip
import sys
from tqdm import tqdm

###### Helper functions ######
def getNumHeaderLines(vcf_filename, num_lines_to_check = 1000):
    num_header_lines = 0
    num_lines_checked = 0
    with gzip.open(vcf_filename, 'rb') as f:
        nextline = f.readline().decode('utf-8')
        num_lines_checked += 1
        while nextline[0:2] == "##" and num_lines_checked <= num_lines_to_check:
            num_header_lines += 1
            nextline = f.readline().decode('utf-8')
            num_lines_checked += 1
    return num_header_lines

def getGT(GT):
    GT1 = 0
    GT2 = 0
    if '/' in GT:
        GT1 = GT.split('/')[0]
        GT2 = GT.split('/')[1]
    elif '|' in GT:
        GT1 = GT.split('|')[0]
        GT2 = GT.split('|')[1]
    try:
        GT1 = int(GT1)
    except:
        GT1 = np.nan
    try:
        GT2 = int(GT2)
    except:
        GT2 = np.nan
    return GT1, GT2

def getAF(GTlist, allele = 1):
    totalGT = 0
    totalALT = 0
    GTnumList = []
    for GT in GTlist:
        GT1, GT2 = getGT(GT)
        GTnumList.append(GT1)
        GTnumList.append(GT2)
    cleanedGT = [x for x in GTnumList if not np.isnan(x)]
    totalGT = len(cleanedGT)
    totalALT = sum([1 for x in cleanedGT if x==allele])
    return float(totalALT / totalGT)

def getDosage(GTlist, allele=1):
    DSlist = []
    for GT in GTlist:
        GT1, GT2 = getGT(GT)
        if np.isnan(GT1) or np.isnan(GT2):
            dosage = np.nan
        else:
            dosage = sum([GT1 == allele, GT2 == allele])
        DSlist.append(dosage)
    return DSlist

def getGP(GTlist, allele=1):
    GPlist = []
    for GT in GTlist:
        GT1, GT2 = getGT(GT)
        if np.isnan(GT1) or np.isnan(GT2):
            GPlist.append('.')
        else:
            dosage = sum([GT1 == allele, GT2 == allele])
            if dosage == 0:
                GPlist.append('1,0,0')
            elif dosage == 1:
                GPlist.append('0,1,0')
            elif dosage == 2:
                GPlist.append('0,0,1')
            else:
                sys.exit(f'Invalid genotype input GT={GT}, dosage={dosage}')
    return GPlist

def appendToGT(GTlist, DSlist, GPlist):
    newGTlist = []
    num_alt_alleles = len(DSlist)
    for GT_index in np.arange(len(GTlist)):
        GT = GTlist[GT_index]
        newGT = GT + ':'
        # append dosage:
        for allele in np.arange(num_alt_alleles):
            DS = DSlist[allele][GT_index]
            if allele > 0:
                newGT += ','
            if np.isnan(DS):
                newGT += '.'
            else:
                newGT +=  str(DS)
        newGT = newGT + ':'
        # append GP:
        for allele in np.arange(num_alt_alleles):
            GP = GPlist[allele][GT_index]
            if allele > 0:
                newGT += ','
            newGT += GP
        newGTlist.append(newGT)
    return newGTlist

def writeHeader(infile, outfilename):
    with gzip.open(infile, 'rb') as f_in:
        with gzip.open(outfilename, 'wb') as f_out:
            nextline = f_in.readline().decode('utf-8')
            while nextline[0:1] == "#":
                f_out.write(bytes(nextline, 'UTF-8'))
                nextline = f_in.readline().decode('utf-8')


##############################
## MAIN
##############################
numHeaderLines = getNumHeaderLines('23-SNPs_to_add_back_normalized.vcf.gz')
vcfchunks = pd.read_csv("23-SNPs_to_add_back_normalized.vcf.gz", compression='gzip', sep='\t', skiprows=numHeaderLines, header=0, chunksize=1000)
#vcfchunk = pd.read_csv("23-SNPs_to_add_back_normalized_chunk.vcf.gz", compression='gzip', sep='\t', skiprows=numHeaderLines, header=0)
outfilename = '24-SNPs_to_add_back_reformatted.vcf.gz'

writeHeader('23-SNPs_to_add_back_normalized.vcf.gz', outfilename)
for vcfchunk in vcfchunks:
    infocol = list(vcfchunk['INFO'])
    formatcol = list(vcfchunk['FORMAT'])
    gt_df = vcfchunk.iloc[:, 9:]
    formatcol = [x.replace('GT','GT:DS:GP') for x in formatcol]
    for row in np.arange(vcfchunk.shape[0]):
        if vcfchunk.iloc[row,3] == '.':
            sys.exit("Missing REF alleles detected! Exiting")
        if vcfchunk.iloc[row,4] == '.':
            infocol[row] = 'DR2=1.00;AF=0'
            gt_df.iloc[row,:] = [str(x) + ':0:1,0,0' for x in gt_df.iloc[row,:]]
            continue
        GTlist = list(gt_df.iloc[row,:])
        alt_alleles = vcfchunk.iloc[row,4].split(',')
        allele_counter = 0
        AF = []
        DSlist = []
        GPlist = []
        for alt_allele in alt_alleles:
            allele_counter += 1            
            AF.append(getAF(GTlist, allele = allele_counter))
            DSlist.append(getDosage(GTlist, allele = allele_counter))
            GPlist.append(getGP(GTlist, allele = allele_counter))
        AFstring = ','.join([f'{x:.4f}' for x in AF])
        infocol[row] = f'DR2=1.00;AF={AFstring}'
        gt_df.iloc[row,:] = appendToGT(GTlist, DSlist, GPlist)
    vcfchunk['INFO'] = infocol
    vcfchunk['FORMAT'] = formatcol
    vcfchunk.iloc[:, 9:] = gt_df
    vcfchunk.to_csv(outfilename, sep='\t', index=False, mode='a', compression='gzip', encoding='utf-8', header=False)

