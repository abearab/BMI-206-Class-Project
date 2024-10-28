

library(ggplot2)
library(dplyr)


make_plot_data = function(data_filename, 
                          plot_chrom, 
                          plot_region_left, 
                         plot_region_right){

if (grepl("uk", data_filename) && !grepl("age_meno", data_filename)) {
  
  chrom_col <- "hm_chrom"
  pval_col <- "p_value"
  pos_col <- "hm_pos"
  rsid_col <- "hm_rsid"
  
} else if (grepl("finn_gen", filename)){
  
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
                  select = c(chrom_col= factor(),
                             pval_col = numeric(),
                             pos_col = integer(),
                             rsid_col = character())
                  # colClasses = c("factor",
                  #                "numeric",
                  #                "integer",
                  #                "character")
                  )

print(head(data))
# rename cols
data = data |> 
       dplyr::rename(pos = !!pos_col,
                     chrom = !!chrom_col,
                     pval = !!pval_col,
                     rsid = !!rsid_col)


print(head(data))

# filtering correct genomic region 
data = data |> 
  dplyr::filter(chrom == plot_chrom) #  #chromosomal region

data = data |> # base regions
  dplyr::mutate(pos = pos / 10^6) |> 
  dplyr::filter(pos >= plot_region_left & pos <= plot_region_right) # convert to mB - and filter proper figure region

# 2. converting p-values for plotting 
data = data |> 
       dplyr::mutate(log_pval = -1*log(pval, base = 10)) 

return(data)
}


# type 1 diabetes is looking correct 

data = make_plot_data(data_filename = "output/gwas_ss_filt/ra_uk_bb.h.filt.tsv",
                      plot_chrom = 6,
                      plot_region_left = 90.1,
                      plot_region_right = 90.4)

data |>
     ggplot(aes(x=pos, y = log_pval)) + 
     geom_point(size = 0.75, alpha = 1) + 
     theme_bw() + 
     labs(x = "Position (chr 6, Mb)", y = "-log_10(P_GWAS)")


