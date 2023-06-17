#script to generate a range +/- 500kb around each transcript to use to assign relevant SNPs
import pandas as pd

#Read in file of transcripts and locations
transcripts = pd.read_csv(snakemake.input[0], sep='\s+', header=0)

#add new columns to transcripts dataframe with range of +/- 500kb
transcripts['range_lower'] = transcripts['s1']-500000

#set any value less than 0 to 0 (i.e. if the gene is within 500kb of the start of the chromosome)
transcripts['range_lower'] = transcripts['range_lower'].clip(lower=0)
transcripts['range_upper'] = transcripts['s2']+500000

#make a directory to store the SNP ranges unless directory already exists
import os
if not os.path.exists(snakemake.output[0]):
    os.mkdir(snakemake.output[0])

#write each row of chr, range_lower / range_upper and geneid to a new file named with geneid:
for index, row in transcripts.iterrows():
    with open(snakemake.output[0] + '/' + row['geneid'] + '_range.txt', 'w') as f:
        f.write(str(row['chr']) + '\t' + str(row['range_lower']) + '\t' + str(row['range_upper']) + '\t' + str(row['geneid']))
