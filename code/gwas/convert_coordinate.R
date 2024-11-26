

hg37_data = c('output/gwas_ss_filt/ra_uk_bb.h.filt.fig4b.tsv',
              'output/gwas_ss_filt/ra_uk_bb.h.filt.fig4a.tsv',
              'output/gwas_ss_filt/hypo_uk_bb.h.filt.fig4a.tsv',
              'output/gwas_ss_filt/endo_uk_bb.h.filt.fig4c.tsv',
              'output/gwas_ss_filt/ovary_cys_uk_bb.h.filt.fig4c.tsv',
              'output/gwas_ss_filt/menorrhagia_uk_bb.h.filt.fig4c.tsv',
              'output/gwas_ss_filt/age_meno_uk_bb.h.filt.fig4c.tsv')

for(data_filename in hg37_data){

df = data.table::fread(data_filename)

df = df |> 
     dplyr::mutate(chrom = paste0("chr", chrom))

# df = df |> 
#      dplyr::mutate(pos = 10^6 * pos)

df = df |> 
      dplyr::mutate(end = pos)

df = df |> 
     dplyr::select(chrom, pos, end, pval)

save_filename = stringr::str_replace(data_filename,
                                     pattern = "gwas_ss_filt",
                                     replacement = "gwas_hg37")

data.table::fwrite(df, 
                   file = save_filename,
                   sep= "\t", 
                   col.names = F
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
