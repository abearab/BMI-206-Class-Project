

module load miniforge3/24.7.1-0
conda activate bmi-206-group
conda install conda-forge::gsutil

# get all the file names for the finemapping data v6

cd data/fine_map && gsutil cp gs://finngen-public-data-r6/summary_stats/R6_manifest.tsv



fingenn_manifest = data.table::fread("data/fine_map/R6_manifest.tsv")
finngenn_cred_filenames = data.table::fread("finngen_filename_cred_sum.tsv")

finngenn_cred_filenames = stringr::remove_string("finngen_R6_", finngenn_cred_filenames)
finngenn_cred_filenames = stringr::remove_string("SUSie")

data/fine_
finngen_R6_C3_

gsutil ls gs://finngen-public-data-r6/finemapping/summaries/ |\
           grep '\.SUSIE\.snp\.filter\.tsv$' |\
           sed 's|.*/||' > data/fine_map/finngen_cred_sum.tsv
