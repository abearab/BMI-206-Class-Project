
{
library(ggplot2)
library(dplyr)
library(LDlinkR)
library(purrr)
}

plot_gg_manhattan = function(data, 
                             data_filename, 
                             top_rsid = "unlabelled",
                             plot_region_left,
                             plot_region_right,
                             breaks,
                             breaks_lab){
  
  message("\n Make ggplot for ", data_filename)
  
  if(top_rsid != "unlabelled"){

  data |>
    ggplot(aes(x=pos, 
               y = log_pval, 
               #fill = R2, 
               color = R2
               )) +
    geom_point(size = 0.5, 
               alpha = 1) +
    xlim(plot_region_left, 
         plot_region_right) + 
      scale_x_continuous(expand = c(0,0),
                         breaks = breaks,
                         labels = breaks_lab) + 
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),   # Remove minor grid lines
      axis.text.x = element_text(size = 7),
      axis.text.y = element_text(size = 7),
      # axis ticks and lines
      axis.ticks = element_line(linewidth = 0.5),
      axis.line.x = element_line(linewidth = 0.25),
      axis.line.y= element_line(linewidth = 0.25),
      # making legend smaller: 
      legend.text = element_text(size = 6),
      legend.title = element_text(size = 7),
      legend.key.size = unit(0.05, "cm")
    ) + 
    labs(x = "Position (chr 6, Mb)", 
         y = expression(-log[10](P[GWAS])) 
    ) + 
    scale_colour_fermenter(name = expression(LD~(r^2)), palette = "YlGnBu") # +
    #scale_fill_fermenter(name = expression(LD~(r^2)), palette = "YlGnBu")
         #y = expression(paste0("-log_10(P_GWAS)")
    
  } else {
    
    data |>
      ggplot(aes(x=pos, y = log_pval)) +
      geom_point(size = 0.5, alpha = 1) +
      xlim(plot_region_left, plot_region_right) + 
      scale_x_continuous(expand = c(0,0),
                         breaks = breaks,
                         labels = breaks_lab) +
      labs(x = "Position (chr 6, Mb)", 
           y = expression(-log[10](P[GWAS])) 
      ) + 
      theme_minimal() +
      theme(
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),   # Remove minor grid lines
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.ticks = element_line(linewidth = 0.5),
        axis.line.x = element_line(linewidth = 0.25),
        axis.line.y= element_line(linewidth = 0.25)
      ) 
     
  }
  
  
  plotfile = stringr::str_replace(data_filename,
                                  pattern = ".tsv",
                                  replacement = ".png")
  
  plotfile = stringr::str_replace(plotfile,
                                  pattern = "output/gwas_plotting_data",
                                  replacement = "figures")
  
  
  message("\n Saving plot as ", plotfile, "\n \n \n")
  
  
  suppressMessages(ggsave(plotfile,
                          width = 20,
                          height = 4,
                          units = "cm"))
}

################ Figure 4a. ############# 
{

fig_4a_data=c(#"output/gwas_plotting_data/ra_uk_bb.h.filt.fig4a.tsv",
             "output/gwas_plotting_data/t1d_uk_bb.h.filt.fig4a.tsv"
              #"output/gwas_plotting_data/hypo_uk_bb.h.filt.fig4a.tsv",
             #"output/gwas_plotting_data/finngen_atopic_derm.filt.fig4a.tsv"
             )

#fig_4a_data = "output/gwas_plotting_data/finngen_atopic_derm.filt.tsv"
# type 1 diabetes is looking correct 


for(data_filename in fig_4a_data){
  
  top_rsid = "rs72928038"
  top_rsid  = "unlabelled"
  breaks = c(90.2, 90.25, 90.3, 90.35)
  breaks_lab = c("90.2", "", "90.3", "")
  plot_region_left = 90.2 - 0.5 * 0.75
  plot_region_right = 90.35 + 0.5 * 0.25

  # data = make_plot_data(data_filename = data_filename,
  #                       plot_chrom = 6,
  #                       plot_region_left = plot_region_left,
  #                       plot_region_right = plot_region_right,
  #                      top_rsid = top_rsid
  #                       )
  # 
  data = data.table::fread(data_filename)
  
  plot_gg_manhattan(data, 
                    data_filename = data_filename,
                    plot_region_left = plot_region_left,
                    plot_region_right = plot_region_right,
                    breaks = breaks,
                    breaks_lab = breaks_lab,
                    top_rsid = top_rsid)

} 

}

################# Figure 4b. ###################### 

{
  
fig_4b_data=c("output/gwas_plotting_data/ra_uk_bb.h.filt.fig4b.tsv",
             "output/gwas_plotting_data/t1d_uk_bb.h.filt.fig4b.tsv"
)

for(data_filename in fig_4b_data){
  
  plot_region_left = 26.05 - 0.05 * 0.75
  plot_region_right = 26.150 + 0.05* 0.25
  breaks = c(26.05, 26.10, 26.15)
  breaks_lab = c("26.05", "26.10", "26.15")
  top_rsid = "rs35944082"
  top_rsid = "unlabelled"
  
  data = data.table::fread(data_filename)
  
  plot_gg_manhattan(data, 
                    data_filename = data_filename,
                    plot_region_left = plot_region_left,
                    plot_region_right = plot_region_right,
                    breaks = breaks,
                    breaks_lab = breaks_lab,
                    top_rsid = top_rsid)
  
} 

}


########### Figure 4c. ############ 

{

fig_4c_data=c("output/gwas_plotting_data/endo_uk_bb.h.filt.fig4c.tsv",
             "output/gwas_plotting_data/ovary_cys_uk_bb.h.filt.fig4c.tsv",
             "output/gwas_plotting_data/menorrhagia_uk_bb.h.filt.fig4c.tsv",
             "output/gwas_plotting_data/age_meno_uk_bb.h.filt.fig4c.tsv"
)

for(data_filename in fig_4c_data){
  
  top_rsid = "rs11031006"
  top_rsid = "unlabelled"
  breaks = c(30.0, 30.1, 30.2, 30.3, 30.4)
  breaks_lab = c("30.0", " ", "30.2", "", "30.4")
  plot_region_left = 30.0 - 0.5 * 0.1 
  plot_region_right = 30.4 + 0.5 * 0.1 
  # data = make_plot_data(data_filename = data_filename,
  #                       plot_chrom = 11,
  #                       plot_region_left = 29.9,
  #                       plot_region_right = 30.4,
  #                       top_rsid = top_rsid)
  # 
  data = data.table::fread(data_filename)
  
  plot_gg_manhattan(data, 
                    data_filename = data_filename,
                    plot_region_left = plot_region_left,
                    plot_region_right = plot_region_right,
                    breaks = breaks,
                    breaks_lab = breaks_lab,
                    top_rsid = top_rsid)
  
} 

}
