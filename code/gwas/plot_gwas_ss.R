

library(ggplot2)
library(dplyr)


make_plot_data = function(data_filename, 
                          plot_chrom, 
                          plot_region_left, 
                         plot_region_right){
  

message("\n Making data for ", data_filename)
  

if (grepl("uk", data_filename) && 
    !grepl("age_meno", data_filename)) {
  
  chrom_col <- "hm_chrom"
  pval_col <- "p_value"
  pos_col <- "hm_pos"
  rsid_col <- "hm_rsid"
  
} else if (grepl("finngen", data_filename)){
  
  chrom_col <- "#chrom"
  pval_col <- "pval"
  pos_col <- "pos"
  rsid_col <- "rsids"
  
} else {
  
  chrom_col <- "chromosome"
  pval_col <- "p_value"
  pos_col <- "base_pair_location"
  rsid_col <- "rsid"
  
}
  
data = data.table::fread(data_filename, 
                  select = list(factor = chrom_col,
                                numeric = pval_col,
                                integer = pos_col,
                                character = rsid_col)
                  )

#print(head(data))

# rename cols
data = data |> 
       dplyr::rename(pos = !!pos_col,
                     chrom = !!chrom_col,
                     pval = !!pval_col,
                     rsid = !!rsid_col)


# filtering correct genomic region 
data = data |> 
  dplyr::filter(chrom == plot_chrom) #  #chromosomal region

print(head(data))

data = data |> # base regions
  dplyr::mutate(pos = pos / 10^6) |> 
  dplyr::filter(pos >= plot_region_left - 1 & 
                pos <= plot_region_right + 1) # convert to mB - and filter proper figure region

#print(head(data))

# 2. converting p-values for plotting 
data = data |> 
       dplyr::mutate(pval = as.numeric(pval))

data = data |> 
       dplyr::mutate(log_pval = -1*log(pval, base = 10)) 

#print(head(data))

return(data)

}


plot_gg_manhattan = function(data){
  
  message("\n Make ggplot for ", data_filename)
  
  data |>
    ggplot(aes(x=pos, y = log_pval)) +
    geom_point(size = 0.5, alpha = 1) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank()   # Remove minor grid lines
    ) + 
    labs(x = "Position (chr 6, Mb)", 
         y = expression(-log[10](P[GWAS])) 
    )
         #y = expression(paste0("-log_10(P_GWAS)")
  
  
  plotfile = stringr::str_replace(data_filename,
                                  pattern = ".filt.tsv",
                                  replacement = ".png")
  
  plotfile = stringr::str_replace(plotfile,
                                  pattern = "output/gwas_ss_filt",
                                  replacement = "figures")
  
  
  message("\n Saving plot as ", plotfile, "\n \n \n")
  
  
  suppressMessages(ggsave(plotfile))
}

################ Figure 4a. ############# 


fig_4a_data=c("output/gwas_ss_filt/ra_uk_bb.h.filt.tsv",
             "output/gwas_ss_filt/t1d_uk_bb.h.filt.tsv",
              "output/gwas_ss_filt/hypo_uk_bb.h.filt.tsv",
             "output/gwas_ss_filt/finngen_atopic_derm.filt.tsv"
           )

#fig_4a_data = "output/gwas_ss_filt/finngen_atopic_derm.filt.tsv"
# type 1 diabetes is looking correct 


for(data_filename in fig_4a_data){

  data = make_plot_data(data_filename = data_filename,
                        plot_chrom = 6,
                        plot_region_left = 90.1,
                        plot_region_right = 90.4)
  
  
  plot_gg_manhattan(data)

} 


################# Figure 4b. ###################### 

fig_4b_data=c("output/gwas_ss_filt/ra_uk_bb.h.filt.tsv",
             "output/gwas_ss_filt/t1d_uk_bb.h.filt.tsv"
)

for(data_filename in fig_4b_data){
  
  data = make_plot_data(data_filename = data_filename,
                        plot_chrom = 4,
                        plot_region_left = 26,
                        plot_region_right = 26.20)
  
  
  plot_gg_manhattan(data)
  
} 


########### Figure 4c. ############ 


fig_4c_data=c("output/gwas_ss_filt/endo_uk_bb.h.filt.tsv",
             "output/gwas_ss_filt/ovary_cys_uk_bb.h.filt.tsv",
             "output/gwas_ss_filt/menorrhagia_uk_bb.h.filt.tsv",
             "output/gwas_ss_filt/age_meno_uk_bb.h.filt.tsv"
)

for(data_filename in fig_4c_data){
  
  data = make_plot_data(data_filename = data_filename,
                        plot_chrom = 11,
                        plot_region_left = 29.9,
                        plot_region_right = 30.4)
  
  
  plot_gg_manhattan(data)
  
} 
