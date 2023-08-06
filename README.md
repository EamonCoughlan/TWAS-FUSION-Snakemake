# TWAS-FUSION-Snakemake
Snakemake workflow for FUSION model generation and TWAS

**These are the files required, with paths mapped in files/twas_config.yaml:**

GWAS summary stats (config gwasfile)

Plink (bed/bim/fam) files of the genotypes from the RNA sequencing, for each condition (naming convention prefix_condition, from config), in a /data subfolder

Plink files of the LD reference, in the /data subfolder (it uses a full genome file, so merge into one if you have split by chr. Workflow will split it later when required.

File with genome location of each gene/transcript from sequencing (4 tab separated columns (ID | chr | startbp | endbp)) (config transcript_coords)

File with expression levels per transcript, mapped data/expression_{condition}.txt




# Notes for users of EDDIE (University of Edinburgh Computing Cluster):

This workflow has hundreds of thousands of small jobs; it is recommended to make a qlogin request for 50+ nodes of 1-2GB each.
It will take around 4 hours to build the DAG (what Snakemake does when initiating/planning workflow)

Recommended run-mode:
`snakemake --cores all --scheduler greedy --keep-going`

If you run without specifying greedy scheduler it defaults to it anyway but adds a 10 or so second delay to each job submission (ie. a lot of time on 750,000 jobs).
--keep-going instructs snakemake to continue if an individual job fails - I think I've now added escapes anyway for failed jobs (eg. a transcript with no variants) to allow continuity in the DAG, but I use keep-going anyway just in case to avoid killing the workflow unnecessarily. Any major issues preventing the run from finishing (eg. missing input files) will still cause it to end.

Before running, run this to avoid filling up /tmp directory on your node:
`source /exports/applications/support/set_qlogin_environment.sh`

It almost certainly won't run in the 48 hours allowed per qlogin session. You can restart the run from wherever it stopped, but use `snakemake --unlock` first or you will get a fail message about a locked directory (after the 4 hour DAG build).

**Documentation - FUSION scripts & packaged software (gcta_nr_robust)**

This workflow uses adapted forms of the FUSION scripts from the gusevlab, with info available here:

https://github.com/gusevlab/fusion_twas

http://gusevlab.org/projects/fusion/
