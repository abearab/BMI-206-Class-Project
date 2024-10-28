

## Steps in GWAS data processing


1. *download_gwas_ss.sh*

Download, and initial chromosome-wide filtering for GWAS full summary statistics. 

Output: 

- data/gwas_ss/* (full summary statistics files)
- output/gwas_ss_filt/* (summary statistics filtered for all variants on chromosome 4,6 & 11)

2. *filter_gwas_ss.sg*

Filtering of GWAS summary statistics for plotting the full GWAS summary statistics. 

Output: 
- output/gwas_ss_filt/* (edited; summary statistics filtered for just the regions and columns needed for plotting)


3. *plot_gwas_ss.R*

Create manhattan plots using the filtered summary statistics produced in step 2. 


4. *download_finemap.sh* 


5. *filter_finemap.sh* 




