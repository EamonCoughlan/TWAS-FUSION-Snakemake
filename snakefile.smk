#Snakemake Workflow to generate TWAS models using FUSION and run TWAS
#Developed by Eamon Coughlan, 2023, using pipeline from Erola Pairo-Castineira.

configfile: 'files/twasfusion_config.yml'
CONDITIONS = config['conditions']
fileprefix = 'data/' + config['seqfile']
#one rule to rule them all (and with Snakemake bind them)
rule all:
    input:
        expand('results/{condition}/TWAS_{i}.dat', condition = CONDITIONS, i = range(1, 23))

#maybe add a rule here to reformat GWAS file to appropriate input format in future

#find SNPs shared between the reference,GWAS and RNAseq data
rule get_shared_SNPs:
    input:
        fileprefix + '_{condition}.bim',
        'data/reference_panel.bim',
        config['gwasfile'],
    output:
        'data/{condition}/shared_SNPs.txt',
    script:
        'scripts/SM-find_common_SNPs.py'

#create merged plink file with only shared SNPs
rule create_merge_bim:
    input:
        fileprefix + '_{condition}.bed',
        fileprefix + '_{condition}.bim',
        fileprefix + '_{condition}.fam',
        'data/{condition}/shared_SNPs.txt'
    output:
        'data/merge_{condition}.bed',
        'data/merge_{condition}.bim',
        'data/merge_{condition}.fam',
    shell:
        'plink --bfile {fileprefix}_{wildcards.condition} --extract data/{wildcards.condition}/shared_SNPs.txt --make-bed --out data/merge_{wildcards.condition}'

#get range +/- 500kb for each transcript and create range file
checkpoint find_SNP_ranges:
    input:
        config['transcript_coords'] + '_{condition}.txt'
    output:
        directory('data/{condition}/snp_ranges/')
    script:
        'scripts/SM-find_SNP_ranges.py'

#get names of range files as snakemake variable {transcriptID} for use in next rule
def get_range_files(wildcards):
    checkpoint_output = checkpoints.find_SNP_ranges.get(**wildcards).output[0]
    return expand('models/{{condition}}/{transcriptID}.wgt.RDat', condition = CONDITIONS, transcriptID = glob_wildcards(os.path.join(checkpoint_output, '{transcriptID}_range.txt')).transcriptID)
    
#create plink files for each gene with the appropriate SNPs from the merge file in the specified range
rule create_beds:
    input:
        'data/{condition}/snp_ranges/{transcriptID}_range.txt',
        'data/merge_{condition}.bed',
        'data/merge_{condition}.bim',
        'data/merge_{condition}.fam',
        'data/{condition}/shared_SNPs.txt'
    output:
        bed = 'data/transcript_beds/{transcriptID}_{condition}.bed',
        bim = 'data/transcript_beds/{transcriptID}_{condition}.bim',
        fam = 'data/transcript_beds/{transcriptID}_{condition}.fam',
    shell:
        'plink --bfile data/merge_{wildcards.condition} --extract range data/{wildcards.condition}/snp_ranges/{wildcards.transcriptID}_range.txt --make-bed --out data/transcript_beds/{wildcards.transcriptID}_{wildcards.condition} || touch {output.bed}; touch {output.bim}; touch {output.fam}'

#merge phenotype info into the fam file
rule incorporate_phenotype_info:
    input:
        'data/transcript_beds/{transcriptID}_{condition}.bed',
        'data/transcript_beds/{transcriptID}_{condition}.bim',
        'data/transcript_beds/{transcriptID}_{condition}.fam',
        'data/expression_{condition}.txt'
    output:
        'data/transcript_beds_with_pheno/{transcriptID}_{condition}.bed',
        'data/transcript_beds_with_pheno/{transcriptID}_{condition}.bim',
        'data/transcript_beds_with_pheno/{transcriptID}_{condition}.fam',
    script:
        'scripts/SM-merge_pheno_info.py'
    #bed and bim files are just copied over to new file with os.system cp - might not work if running on Windows

#make directories models/* and tmp_*/ for each condition
import os
for condition in CONDITIONS:
    if not os.path.exists('tmp_' + condition):
        os.mkdir('tmp_' + condition)
    if not os.path.exists('models/'):
        os.mkdir('models/')
    if not os.path.exists('models/' + condition):
        os.mkdir('models/' + condition)

#generate models using FUSION
rule calculate_models:
    input:
        bed = 'data/transcript_beds_with_pheno/{transcriptID}_{condition}.bed',
        bim = 'data/transcript_beds_with_pheno/{transcriptID}_{condition}.bim',
        fam = 'data/transcript_beds_with_pheno/{transcriptID}_{condition}.fam',
    output:
        hsq = 'models/{condition}/{transcriptID}.hsq',
        RDat = 'models/{condition}/{transcriptID}.wgt.RDat'
    shell:
        'Rscript scripts/FUSION.compute_weights.R --bfile data/transcript_beds_with_pheno/{wildcards.transcriptID}_{wildcards.condition} --tmp tmp_{wildcards.condition}/ --hsq_p 1 --verbose 1 --save_hsq --models lasso,top1,enet --out models/{wildcards.condition}/{wildcards.transcriptID} || touch {output.hsq}; touch {output.RDat}'

#filler rule to create files which were not completed by earlier rules (incorporated into previous to avoid out-competing the pipeline that actually generates files)
#rule filler:
    #input:
        #'data/{condition}/snp_ranges/{transcriptID}_range.txt'
    #output:
        #'models/{condition}/{transcriptID}.wgt.RDat'
    #shell:
        #'touch {output}'

rule collate_models:
    input:
        get_range_files
    output:
        'data/{condition}/model_list.txt'
    shell:
        r'''for file in {input}; do
                if [ -s ${{file}} ]; then
                    echo ${{file}}
                    echo ${{file}} >> {output}
                else
                    echo ${{file}} is empty
                    continue
                fi
            done'''

#create profile (summary of model characteristics) for each gene
rule create_profile:
    input:
        'data/{condition}/model_list.txt'
    output:
        'data/{condition}.profile'
    script:
        'scripts/SM-createprofile.R'

#create pos file listing which models to use in the TWAS
rule create_pos_file:
    input:
        'data/{condition}.profile',
        config['transcript_coords'] + '_{condition}.txt',
        fileprefix + '_{condition}.fam',
    params:
        pval = config['pvalue']
    output:
        'data/{condition}_' + str(config['pvalue']) + '.pos'
    script:
        'scripts/SM-createpos.py'

#make directory for refpanel chromosomes if not already present
if not os.path.exists('data/reference_panel_chr/'):
    os.mkdir('data/reference_panel_chr/')

#rule to make the reference panel into separate chromosomes, then use those as input to the TWAS rule instead of that chr.bim
rule split_refpanel_by_chromosome:
    input:
       'data/reference_panel.bed',
       'data/reference_panel.bim',
       'data/reference_panel.fam',
    output: 
        'data/reference_panel_chr/{i}.bed',
        'data/reference_panel_chr/{i}.bim',
        'data/reference_panel_chr/{i}.fam',
    shell:
        'plink --bfile data/reference_panel --chr {wildcards.i} --make-bed --out data/reference_panel_chr/{wildcards.i}'

#run TWAS
rule run_TWAS:
    input:
        'data/reference_panel_chr/{i}.bed',
        'data/reference_panel_chr/{i}.bim',
        'data/reference_panel_chr/{i}.fam',
        gwas = config['gwasfile'],
        posfile = 'data/{condition}_' + str(config['pvalue']) + '.pos',
    params:
        pval = config['pvalue']
    output:
        'results/{condition}/TWAS_{i}.dat'
    shell:
        'Rscript scripts/FUSION.assoc_test.R --sumstats {input.gwas} --weights {input.posfile} --weights_dir ./ --ref_ld_chr data/reference_panel_chr/ --chr {wildcards.i} --out {output}'