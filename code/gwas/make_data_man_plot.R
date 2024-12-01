
{
  suppressMessages(library(dplyr))
  library(LDlinkR)
  library(purrr)
}

# making plot data takes ~5 - 10 mins 
# due to LD labeling 
make_plot_data = function(data_filename, 
                          plot_chrom, 
                          plot_region_left, 
                          plot_region_right, 
                          top_rsid = "unlabelled",
                          top_rsid_pos){
  
  
  message("\n Making data for ", data_filename)
  
  
  if (grepl("uk", data_filename) && 
      grepl("t1d", data_filename)) {
    
    chrom_col <- "hm_chrom"
    pval_col <- "p_value"
    pos_col <- "hm_pos"
    rsid_col <- "hm_rsid"
    
  } else if (grepl("finngen", data_filename)){
    
    chrom_col <- "#chrom"
    pval_col <- "pval"
    pos_col <- "pos"
    rsid_col <- "rsids"
    
  } else if (grepl("uk", data_filename)){
    
    rsid_col <- NULL
    pval_col <- "pval"
    chrom_col <- "chrom"
    pos_col <- "pos"
    
  } else {
    
    chrom_col <- "chromosome"
    pval_col <- "p_value"
    pos_col <- "base_pair_location"
    rsid_col <- "rsid"
  }
  
  
  
  if (!grepl("gwas_plotting_data", data_filename)){
  
  data = data.table::fread(data_filename, 
                           select = list(
                                         character = rsid_col,            
                                         factor = chrom_col,
                                         integer = pos_col,
                                         numeric = pval_col
                                         )
  )
  
  if(is.null(rsid_col)){
    
  data = data |> 
         tidyr::unite("rsid", 
                      chrom_col:pos_col, 
                      remove = F, 
                      sep = ":",
                      na.rm = T)
  
  rsid_col = "rsid"
  
  if(top_rsid != "unlabelled"){
    
    top_rsid = top_rsid_pos
    
  }
  }
  
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
  
  
  data = data |> # base regions
    dplyr::mutate(pos = pos / 10^6) |> 
    dplyr::filter(pos >= plot_region_left - 0 & 
                    pos <= plot_region_right + 0) # convert to mB - and filter proper figure region
  
  #print(head(data))
  
  # 2. converting p-values for plotting 
  data = data |> 
    dplyr::mutate(pval = as.numeric(pval))
  
  data = data |> 
    dplyr::mutate(log_pval = -1*log(pval, base = 10)) 
  
  } else {
    
    data = data.table::fread(data_filename)
  }
  # 
  print(head(data))
  
  # get relevant LD info
  
  if(top_rsid != "unlabelled"){
    
    message("Getting LD data for ", top_rsid)
    
    if(grepl(":", data[1,] |> dplyr::pull("rsid"))){
      
      top_rsid_keep = top_rsid
      top_rsid = top_rsid_pos
    }
    
    top_snp_pos = data |> 
      dplyr::filter(rsid == top_rsid) |> 
      dplyr::pull("pos")
    
    if(length(top_snp_pos) == 0){
      message("snp not in dataset")
    }
    
    if(length(top_snp_pos) != 0){
      
      if(grepl(":",top_rsid)){
        data = data |> 
               dplyr::mutate(rsid = paste0("chr", rsid))
        
        top_rsid = paste0("chr", top_rsid)
        
        genome_build = "grch38"
      } else {
        genome_build = "grch37"
      }
    
    #message("\n for ", length(top_snp_pos), " snps")
    
    win_size = max(abs(plot_region_right - top_snp_pos),
                   abs(plot_region_left  - top_snp_pos))
    
    print(win_size)
    
    
    high_ld = suppressMessages(
      LDlinkR::LDproxy(snp = top_rsid, 
                       pop = "EUR", 
                       r2d = "r2", 
                       genome_build = genome_build, #"grch37",
                       win_size = win_size*10^6,
                       token = Sys.getenv("LDLINK_TOKEN")
      )
    )
    
    print(names(high_ld))
    print(head(high_ld))
    
    if(!grepl(":",top_rsid)){
    high_ld = high_ld |> 
      dplyr::select(dplyr::all_of("RS_Number", "R2")) |> 
      dplyr::rename(rsid = RS_Number)
    } else {
      high_ld = high_ld |> 
        dplyr::select(dplyr::all_of(c("Coord", "R2"))) |> 
        dplyr::rename(rsid = Coord)
    }
    
    low_ld = setdiff(data$rsid, high_ld$rsid) |> 
      unique()
    
    # filtering out non labelled rsids
    low_ld = low_ld[low_ld != " . "]
    low_ld = low_ld[!grepl(pattern = "\\.", low_ld)]
    if(!grepl(":",top_rsid)){
    low_ld = low_ld[grepl(pattern = "rs", low_ld)]
    
    }
    low_ld_split <- split(low_ld, 
                          ceiling(seq_along(low_ld) / 200)
    )
    
    # mat = suppressMessages(LDlinkR::LDmatrix(
    #   snps = c(chunk, top_rsid),
    #   pop = "EUR",
    #   r2d = "r2",
    #   genome_build = "grch37",
    #   token = Sys.getenv("LDLINK_TOKEN")
    # )
    
    if(exists("top_snp_pos")){
      
      browser()
    
    low_ld_mat <- purrr::map_dfr(low_ld_split[1:2], function(chunk) {
      mat = suppressMessages(LDlinkR::LDmatrix(
        snps = c(chunk, top_rsid),
        pop = "EUR",
        r2d = "r2",
        genome_build = genome_build,
        token = Sys.getenv("LDLINK_TOKEN")
      )[,c("RS_number",  top_rsid_keep)]
      )
    },
    .progress = T
    )
    
    = SNPlocs.Hsapiens.dbSNP150.GRCh38::rsid2loc(rsids, caching=TRUE)
    
    
    rsid_convert = low_ld_mat |> 
                   dplyr::pull("RS_number") |> 
                   unique()
    
    rsid_convert = split(rsid_convert, 
                         ceiling(seq_along(rsid_convert) / 5000
                         )
                                 )
    
    
    
    rsid_mapping <- purrr::map_dfr(rsid_convert , function(chunk) {
      suppressMessages(LDlinkR::SNPchip(
        snps = c(chunk),
   #     pop = "EUR",
        chip = "A_SNP5.0",
        genome_build = genome_build,
        token = Sys.getenv("LDLINK_TOKEN")
      )[,c("RS_number",  Position_GRCh38)]
      )
    },
    .progress = T
    )
    
    low_ld_mat = dplyr::inner_join(low_ld_mat, 
                                   rsid_mapping, 
                                   by = "RS_number")
    
    low_ld_mat = low_ld_mat |>
                 dplyr::select(-RS_number)
    
    low_ld_mat = low_ld_mat |> 
                 dplyr::rename(RS_number = Position_GRCh38)
    
    } else {
    
    low_ld_mat <- purrr::map_dfr(low_ld_split, function(chunk) {
      mat = suppressMessages(LDlinkR::LDmatrix(
        snps = c(chunk, top_rsid),
        pop = "EUR",
        r2d = "r2",
        genome_build = genome_build,
        token = Sys.getenv("LDLINK_TOKEN")
      )[,c("RS_number",  top_rsid)]
      )
    },
    .progress = T
    )
    
    }
    
    #if(!grepl(":",top_rsid)){
      
    low_ld = low_ld_mat |>         
      dplyr::rename(rsid = RS_number,
                    R2  = !!top_rsid)
    
    #}
    
    
    
    all_ld = dplyr::bind_rows(high_ld, low_ld) |>
      distinct()
    
    
    
    # low_ld <- split(low_ld, 
    #                 cut(seq_along(low_ld), 
    #                 4, labels = FALSE)
    #                 )
    # 
    # safe_LDpair <- purrr::possibly(function(snp) {
    # 
    # 
    #   suppressMessages(
    #                    LDpair(top_rsid,
    #                    snp,
    #                   pop = "EUR",
    #                    token = Sys.getenv("LDLINK_TOKEN"),
    #                                                      genome_build = "grch38" )
    #                   )[, "r2"] },
    #                                otherwise = NA)
    
    # tic("LDpair")
    # 
    # # low_ld_snps = low_ld
    # # 
    # # low_ld <- purrr:::map_vec(low_ld_snps,
    # #                       safe_LDpair,
    # #                       .progress = T)
    # # 
    # # low_ld = data.frame(rsid = low_ld_snps,
    # #                     R2 = low_ld)
    # 
    # toc()
    # stop()
    

    data = dplyr::left_join(data,
                            all_ld,
                            by = "rsid")
    
    if(grepl("chr",top_rsid)){
    
    data = data |> 
           mutate(rsid = stringr::str_remove(rsid, "chr"))
    }
    
    }
  }
  
  print(head(data))
  
  return(data)
  
}


########### Figure 4a. ############ 

fig_4a_data=c("output/gwas_hg38_preplotting/ra_uk_bb.h.filt.fig4a.tsv",
 "output/gwas_ss_filt/t1d_uk_bb.h.filt.fig4a.tsv",
  "output/gwas_hg38_preplotting/hypo_uk_bb.h.filt.fig4a.tsv",
"output/gwas_ss_filt/finngen_atopic_derm.filt.fig4a.tsv"
)


for(data_filename in fig_4a_data){
  
  top_rsid = "rs72928038"
  top_rsid_pos = "6:90267049"
  top_rsid = "unlabelled"
  plot_region_left = 90.2 - 0.5 * 0.75
  plot_region_right = 90.35 + 0.5 * 0.25
  
  data = make_plot_data(data_filename = data_filename,
                        plot_chrom = 6,
                        plot_region_left = plot_region_left,
                        plot_region_right = plot_region_right,
                        top_rsid = top_rsid,
                        top_rsid_pos = top_rsid_pos
  )
  
  save_filename = stringr::str_replace(data_filename,
                                       pattern = "gwas_ss_filt",
                                       replacement = "gwas_plotting_data")
  
  
  save_filename = stringr::str_replace(data_filename,
                                       pattern = "gwas_hg38_preplotting",
                                       replacement = "gwas_plotting_data")  
  data.table::fwrite(data, 
                     file = save_filename)
  
  
  
} 

################# Figure 4b. ###################### 

{
  
  fig_4b_data=c("output/gwas_hg38_preplotting/ra_uk_bb.h.filt.fig4b.tsv",
                "output/gwas_ss_filt/t1d_uk_bb.h.filt.fig4b.tsv"
  )
  
  for(data_filename in fig_4b_data){
    
    
    top_rsid = "rs35944082"
    top_rsid_pos = "4:26093692"
    top_rsid = "unlabelled"
    plot_region_left = 26.05 - 0.05 * 0.75
    plot_region_right = 26.150 + 0.05* 0.25
    
    data = make_plot_data(data_filename = data_filename,
                          plot_chrom = 4,
                          plot_region_left = plot_region_left,
                          plot_region_right = plot_region_right,
                          top_rsid = top_rsid, 
                          top_rsid_pos = top_rsid_pos)
    
  save_filename = stringr::str_replace(data_filename,
                         pattern = "gwas_ss_filt",
                         replacement = "gwas_plotting_data")
  
  save_filename = stringr::str_replace(data_filename,
                                       pattern = "gwas_hg38_preplotting",
                                       replacement = "gwas_plotting_data")
  
  data.table::fwrite(data, 
                     file = save_filename)
    
    
  } 
  
}


########### Figure 4c. ############ 

  {
  
  fig_4c_data=c("output/gwas_hg38_preplotting/endo_uk_bb.h.filt.fig4c.tsv",
                "output/gwas_hg38_preplotting/ovary_cys_uk_bb.h.filt.fig4c.tsv",
                "output/gwas_hg38_preplotting/menorrhagia_uk_bb.h.filt.fig4c.tsv",
                "output/gwas_hg38_preplotting/age_meno_uk_bb.h.filt.fig4c.tsv"
  )
  
  for(data_filename in fig_4c_data){
    
    top_rsid = "rs11031006"
    top_rsid_pos = "11:30204981"
    top_rsid = "unlabelled"
    plot_region_left = 30.0 - 0.5 * 0.1 
    plot_region_right = 30.4 + 0.5 * 0.1 
    
    data = make_plot_data(data_filename = data_filename,
                          plot_chrom = 11,
                          plot_region_left = plot_region_left,
                          plot_region_right = plot_region_right,
                          top_rsid = top_rsid,
                          top_rsid_pos = top_rsid_pos)

    save_filename = stringr::str_replace(data_filename,
                                         pattern = "gwas_hg38_preplotting",
                                         replacement = "gwas_plotting_data")
    
    
    data.table::fwrite(data, 
                       file = save_filename)
    
  } 
  
}


########### Figure 4a. LD!!!! ############ 

fig_4a_data=c("output/gwas_plotting_data/ra_uk_bb.h.filt.fig4a.tsv",
              "output/gwas_plotting_data/t1d_uk_bb.h.filt.fig4a.tsv",
              "output/gwas_plotting_data/hypo_uk_bb.h.filt.fig4a.tsv",
              "output/gwas_plotting_data/finngen_atopic_derm.filt.fig4a.tsv"
)


for(data_filename in fig_4a_data){
  
  top_rsid = "rs72928038"
  top_rsid_pos = "6:90267049"
  top_rsid = "unlabelled"
  plot_region_left = 90.2 - 0.5 * 0.75
  plot_region_right = 90.35 + 0.5 * 0.25
  
  data = make_plot_data(data_filename = data_filename,
                        plot_chrom = 6,
                        plot_region_left = plot_region_left,
                        plot_region_right = plot_region_right,
                        top_rsid = top_rsid,
                        top_rsid_pos = top_rsid_pos
  )
  
  save_filename = stringr::str_replace(data_filename,
                                       pattern = "gwas_plotting_data",
                                       replacement = "gwas_plotting_data/ld")
  
 
  data.table::fwrite(data, 
                     file = save_filename)
  
  
} 

################# Figure 4b. LD !!! ###################### 

{
  
  fig_4b_data=c("output/gwas_hg38_preplotting/ra_uk_bb.h.filt.fig4b.tsv",
                "output/gwas_ss_filt/t1d_uk_bb.h.filt.fig4b.tsv"
  )
  
  for(data_filename in fig_4b_data){
    
    
    top_rsid = "rs35944082"
    top_rsid_pos = "4:26093692"
    top_rsid = "unlabelled"
    plot_region_left = 26.05 - 0.05 * 0.75
    plot_region_right = 26.150 + 0.05* 0.25
    
    data = make_plot_data(data_filename = data_filename,
                          plot_chrom = 4,
                          plot_region_left = plot_region_left,
                          plot_region_right = plot_region_right,
                          top_rsid = top_rsid, 
                          top_rsid_pos = top_rsid_pos)
    
    save_filename = stringr::str_replace(data_filename,
                                         pattern = "gwas_plotting_data",
                                         replacement = "gwas_plotting_data/ld")
    
    
    data.table::fwrite(data, 
                       file = save_filename)
    
    
  } 
  
}


########### Figure 4c. LD!!! ############ 

{
  
  fig_4c_data=c("output/gwas_plotting_data/endo_uk_bb.h.filt.fig4c.tsv",
                "output/gwas_plotting_data/ovary_cys_uk_bb.h.filt.fig4c.tsv",
                "output/gwas_plotting_data/menorrhagia_uk_bb.h.filt.fig4c.tsv",
                "output/gwas_plotting_data/age_meno_uk_bb.h.filt.fig4c.tsv"
  )
  
  for(data_filename in fig_4c_data){
    
    top_rsid = "rs11031006"
    #top_rsid = "11:30204981"
    top_rsid_pos = "11:30204981"
    #top_rsid = "unlabelled"
    plot_region_left = 30.0 - 0.5 * 0.1 
    plot_region_right = 30.4 + 0.5 * 0.1 
    
    data = make_plot_data(data_filename = data_filename,
                          plot_chrom = 11,
                          plot_region_left = plot_region_left,
                          plot_region_right = plot_region_right,
                          top_rsid = top_rsid,
                          top_rsid_pos = top_rsid_pos)
    
    save_filename = stringr::str_replace(data_filename,
                                         pattern = "gwas_plotting_data",
                                         replacement = "gwas_plotting_data/ld")
    
    
    data.table::fwrite(data, 
                       file = save_filename)
    
  } 
  
}

