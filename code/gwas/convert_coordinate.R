hg37_data = c('output/gwas_ss_filt/ra_uk_bb.h.filt.fig4b.tsv',
              'output/gwas_ss_filt/ra_uk_bb.h.filt.fig4a.tsv',
              'output/gwas_ss_filt/hypo_uk_bb.h.filt.fig4a.tsv',
              'output/gwas_ss_filt/endo_uk_bb.h.filt.fig4c.tsv',
              'output/gwas_ss_filt/ovary_cys_uk_bb.h.filt.fig4c.tsv',
              'output/gwas_ss_filt/menorrhagia_uk_bb.h.filt.fig4c.tsv',
              'output/gwas_ss_filt/age_meno_uk_bb.h.filt.fig4c.tsv')

# # Ensure output directory exists
# if (!dir.exists("output/gwas_hg37")) {
#   dir.create("output/gwas_hg37", recursive = TRUE)
# }

for (data_filename in hg37_data) {
  print(data_filename)
  
  # Read the data
  df = data.table::fread(data_filename)
  
  #Transform the data
  df = df |>
    dplyr::mutate(chrom = paste0("chr", chrom),
                  end = pos) |>
    dplyr::select(chrom, pos, end, pval)
  
  df = df |> 
       dplyr::filter(!is.na(pval))

  # Generate the save filename
  save_filename = stringr::str_replace(
    data_filename,
    pattern = "gwas_ss_filt",
    replacement = "gwas_hg37"
  )

  print("saving as ...")
  print(save_filename)

  # Write the transformed data
  data.table::fwrite(
    df,
    file = save_filename,
    sep = "\t",
    col.names = FALSE
  )
}


# use the ucsc liftover

for(data_filename in hg37_data){

converted = stringr::str_replace(data_filename,
                                     pattern = "gwas_ss_filt",
                                     replacement = "gwas_hg38_converted")

converted = stringr::str_replace(converted,
                                 pattern = ".tsv",
                                 replacement = ".bed")

converted = data.table::fread(converted)

names(converted) = c("chrom", "pos", "end", "pval")

converted = converted |> 
            dplyr::select(-end)

converted = converted |> 
            dplyr::mutate(chrom = stringr::str_remove(chrom,
                                                      pattern = "chr"))

save_filename = stringr::str_replace(data_filename,
                                     pattern = "gwas_ss_filt",
                                     replacement = "gwas_hg38_preplotting")


data.table::fwrite(converted, 
                   file = save_filename,
                   sep= "\t"
)
}
