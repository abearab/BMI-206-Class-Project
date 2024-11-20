

## Steps in GWAS data processing for Manhattan Plots 


<h2> 1. *download_gwas_ss.sh* </h2>

Download, and initial chromosome-wide filtering for GWAS full summary statistics. 

Output: 
- data/gwas_ss/* (full summary statistics files)
- output/gwas_ss_filt/* (summary statistics filtered for all variants on chromosome 4,6 & 11)

<h2> 2. *filter_gwas_ss.sg* </h2>

Further filter the full GWAS summary statistics for the specific regions that are plotted in the Manhattan plots. 

Input:
- output/gwas_ss_filt/*

Output: 
- output/gwas_ss_filt/* (edited; summary statistics filtered for just the regions and columns needed for plotting)

<h2> 3. *make_data_man_plot.sh* (calls *make_data_man_plot.R*) </h2>

Reform the filtered GWAS summary statistics, and integrate them with LD information (from 1000 genomes).
This creates the data.frame that can then directly be plotted. 

Input:
- output/gwas_ss_filt/*

Output: 
- output/gwas_plotting_data/*


<h2> 4. *plot_gwas_ss.R* </h2>

Create Manhattan plots using the datasets created in prior step. 

Input:
- output/gwas_plotting_data/*

Output: 
- figures/*


