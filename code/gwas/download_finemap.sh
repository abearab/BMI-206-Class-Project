

module load miniforge3/24.7.1-0
conda activate bmi-206-group
#conda install conda-forge::gsutil

# get all the file names for the finemapping data v6

gsutil cp gs://finngen-public-data-r6/finemapping/summaries/finngen_R6_ATOPIC_STRICT.SUSIE.snp.filter.tsv "test.v2.tsv"

#gs://finngen-public-data-r6/finemapping/full/finngen_R6_AB1_AMOEBIASIS.SUSIE.snp.bgz

# so instead - let's extract those nanmes
# for full SUSIE results 
# gsutil ls gs://finngen-public-data-r6/finemapping/full |\
#            grep '\.SUSIE\.snp\.bgz$' |\
#            sed 's|.*/||' > data/fine_map/finngen_full/finngen_full_list.tsv

# for snp-filter 
# snps that are part of good quality credible sets as reported in PHENONAME.cred.summary.tsv file.
gsutil ls gs://finngen-public-data-r6/finemapping/summaries/ |\
           grep '\.SUSIE\.snp\.filter\.tsv$' |\
           sed 's|.*/||' > data/fine_map/finngen_snp_filter_list.tsv



# these file names don't match what's avaliable in the bucket
gsutil cp gs://finngen-public-data-r6/summary_stats/R6_manifest.tsv "data/fine_map/R6_manifest.tsv"
        
           
RScript -e '

R6_all = data.table::fread("data/fine_map/R6_manifest.tsv");

traits = R6_all |> 
dplyr::filter(n_cases > 1000) |>
dplyr::pull("phenocode");

#traits = data.table::fread("data/fine_map/R6_manifest.tsv")$phenocode;
filenames = readLines("data/fine_map/finngen_snp_filter_list.tsv");

matching_files <- filenames[sapply(filenames, function(file) {
  any(sapply(traits, function(trait) grepl(trait, file)))
})]

matching_files = data.frame(matching_files)

data.table::fwrite(matching_files,
                    "data/fine_map/finngen_snp_filter_list.tsv",
                    col.names=F
                  )
'

#module load htslib/1.21

# download + basic process all the relevant files
while read -r filename; do

  if [ ! -f "data/fine_map/$filename.pip0.5.nonsyn" ]; then
      # Download each file from Google Cloud Storage and save it to data/fine_map
      gsutil cp gs://finngen-public-data-r6/finemapping/summaries/$filename "data/fine_map/$filename"
  
  # # unzip file
  # bgzip -d "data/fine_map/$filename"
  # filename="${filename%.bgz}"
  
  # filter out snps with PIP <= 0.5
      awk 'NR == 1 {print; next} $5 > 0.5' "data/fine_map/$filename" > "data/fine_map/$filename.pip0.5"

      awk '!($14 ~ /missense_variant|synonymous_variant|frameshift_variant|splice_acceptor_variant|splice_donor_variant|stop_gained/)' "data/fine_map/$filename.pip0.5" > "data/fine_map/$filename.pip0.5.nonsyn"

  #awk '$16 != -1' "data/fine_map/finngen_full/$filename" > "data/fine_map/finngen_full/$filename.cs"
  fi
  
rm "data/fine_map/$filename"
rm "data/fine_map/$filename.pip0.5"
  
  done < data/fine_map/finngen_snp_filter_list.tsv

# combine files if they contain at least one snp 


# Create an empty list to hold filenames
output_file="data/fine_map/all_finngen_trait_cs_filtered.tsv"

> "$output_file"  # Initialize the output file (clear it if it exists)

# Loop through each .snp.filter.tsv file in the folder

for file in data/fine_map/*.pip0.5.nonsyn; do
  # Check if the file has more than 1 line
  line_count=$(wc -l < "$file")
  
  if [ "$line_count" -gt 1 ]; then
      cat "$file" >> "$output"
      echo "" >> "$output"  # Add a newline after each file for separation
  fi
done

# Filtering: step 1

#  selected any traits with case count >1,000


# Filtering: step 2: 

# noncoding fine-mapped loci that did not include any nonsynonymous or splicing variants 

# missense_variant
# synonymous_variant
# frameshift_variant
# splice_acceptor_variant
# splice_donor_variant
# stop_gained

# ? non_coding_transcript_exon_variant
# ? intergenic_variant


# fingenn_manifest = data.table::fread("data/fine_map/R6_manifest.tsv")
# finngenn_cred_filenames = data.table::fread("finngen_filename_cred_sum.tsv")
# 
# finngenn_cred_filenames = stringr::remove_string("finngen_R6_", finngenn_cred_filenames)
# finngenn_cred_filenames = stringr::remove_string("SUSie")

# data/fine_
# finngen_R6_C3_

# gsutil ls gs://finngen-public-data-r6/finemapping/summaries/ |\
#            grep '\.SUSIE\.snp\.filter\.tsv$' |\
#            sed 's|.*/||' > data/fine_map/finngen_cred_sum.tsv
