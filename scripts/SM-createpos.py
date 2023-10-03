import pandas as pd

#import .profile file
profile = pd.read_csv(snakemake.input[0], header=0, delim_whitespace=True)

#import location file
location = pd.read_csv(snakemake.input[1], header=0, delim_whitespace=True)

#import fam file
fam = pd.read_csv(snakemake.input[2], header=None, delim_whitespace=True)
n = len(fam) #number of samples

#filter .profile file for values with pv > 0.05
profile = profile[profile['hsq.pv'] < snakemake.params['pval']]

#extract ENS id from between out_UT/ and .snps in profile id column - not necessary with current format of .profile
#profile['ensid'] = profile['id'].str.extract('out_UT/(.*).snps', expand=True)

#create pos table from location file for all geneid with pv > 0.05 in the .profile file
#loc = location[location['geneid'].isin(profile['ensid'])]
loc = location[location['geneid'].isin(profile['id'])]
#print(loc)

profile['PANEL'] = snakemake.wildcards['condition']
#merge pos with profile based on ensid/geneid
profile = pd.merge(profile, loc, left_on='id', right_on='geneid')

#make dataframe column consisting of id column with 'models/' and {condition} prefix added
profile['path'] = 'models/' + snakemake.wildcards['condition'] + '/' + profile['id'] + '.wgt.RDat'
#make outfile dataframe with columns from profile (PANEL, id, ensid, chr, s1, s2)
#outfile = profile[['PANEL', 'id', 'ensid', 'chr', 's1', 's2']]
outfile = profile[['PANEL', 'path', 'id', 'chr', 's1', 's2']]
outfile['N'] = str(n)

#rename columns to PANEL WGT ID CHR P0 P1 N
outfile.columns = ['PANEL', 'WGT', 'ID', 'CHR', 'P0', 'P1', 'N']

#write outfile to snakemake output
outfile.to_csv(snakemake.output[0], sep='\t', index=False)
