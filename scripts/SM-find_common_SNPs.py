#Script to find SNPs in common across 3 input bim files
import pandas as pd

#read in the 3 files
file1 = pd.read_csv(snakemake.input[0], sep='\s+', header=None) #imputed genotypes from RNAseq
file2 = pd.read_csv(snakemake.input[1], sep='\s+', header=None) #reference genomes
file3 = pd.read_csv(snakemake.input[2], sep='\s+', header=0) #gwas metadata

#create a new dataframe consisting of the intersection of column 2 (SNP ID) from the 3 files
merge12 = pd.merge(file1, file2, on=1)

#merge the new dataframe with 3rd file based on column 2 of the new dataframe and column 1 ('SNP') of 3rd file (GWAS metadata)
common_SNPs = pd.merge(merge12, file3, left_on=1, right_on='SNP')

#write column 2 of common_SNPs to a new file
common_SNPs[1].to_csv(snakemake.output[0], sep='\t', header=False, index=False)
