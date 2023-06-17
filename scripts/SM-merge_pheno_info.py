#import files (bed/bim/fam per transcript and expression file)
import pandas as pd

expr = pd.read_csv(snakemake.input[3], sep='\s+', header=0)

import os
if not os.path.exists('data/transcript_beds_with_pheno/'):
    os.mkdir('data/transcript_beds_with_pheno/')

#get row corresponding to transcriptID from expression file
expr_row = expr.loc[expr['id'] == str(snakemake.wildcards['transcriptID'])]

#add expression values from expr_row to each row of column 6 in fam (if fam file is empty, copy over as with bim and bed files)
try:
    fam = pd.read_csv(snakemake.input[2], sep='\s+', header=None)
    for i in range (0, len(fam[5])):
        fam[5][i] = expr_row.iloc[0,1:][i]
    #write fam to new file
    fam.to_csv(snakemake.output[2], sep=' ', header=False, index=False)
except:
    print('an error occurred, probably due to empty (placeholder) fam file. Safe to ignore')
    os.system('cp ' + snakemake.input[2] + ' ' + snakemake.output[2])



#copy bed and bim files to new output files with os.system
os.system('cp ' + snakemake.input[0] + ' ' + snakemake.output[0])
os.system('cp ' + snakemake.input[1] + ' ' + snakemake.output[1])