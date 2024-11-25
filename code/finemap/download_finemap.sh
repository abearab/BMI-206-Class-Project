#conda activate bmi-206-group
mkdir -p output/fine_map

########### Part 1: Determine: what finngen files to dowload ##################

echo "Determine what FinnGen Files to Download"
echo " "

# a. Get all the file names for the finemapping data v6

# get a full list of avaliable snp-filter trait files
# snp-filter susie datasets are 
# snps that are part of good quality credible sets as reported in PHENONAME.cred.summary.tsv file.

if [ ! -f "data/fine_map/finngen_snp_filter_list.tsv" ]; then
   
   gsutil ls gs://finngen-public-data-r6/finemapping/summaries/ |\
           grep '\.SUSIE\.snp\.filter\.tsv$' |\
           sed 's|.*/||' > data/fine_map/finngen_snp_filter_list.tsv
fi

"data/fine_map/finngen_snp_filter_list.tsv"

# b. the traits tested/file names don't exactly match what's avaliable in the bucket for the manifest file- 
# so let's overlap them
if [ ! -f "data/fine_map/R6_manifest.tsv" ]; then
   gsutil cp gs://finngen-public-data-r6/summary_stats/R6_manifest.tsv "data/fine_map/R6_manifest.tsv"
fi
        
  
# we overlap them while doing filtering step 1: filter for traits with cases > 1,000         
Rscript -e '
Rscript -e '

R6_all = data.table::fread("data/fine_map/R6_manifest.tsv");

traits = R6_all |> 
dplyr::filter(n_cases > 1000) |>
dplyr::pull("phenocode");

#traits = data.table::fread("data/fine_map/R6_manifest.tsv")$phenocode;
filenames = readLines("data/fine_map/finngen_snp_filter_list.tsv");

matching_files <- filenames[sapply(filenames, function(file) {
  any(sapply(traits, function(trait) grepl(trait, file)))
})]

matching_files = data.frame(matching_files);

print(nrows(matching_files))

data.table::fwrite(matching_files,
                   here::here("data/fine_map/finngen_snp_filter_list.tsv"),
                    col.names=F
                  )
'

exit

########### Step 2. Download the files, and do basic processing ##############

#module load htslib/1.21

echo "Download + do basic filtering of the relevant files"
echo " "
echo " "

# download + basic process all the relevant files
while read -r filename; do

  if [ ! -f "data/fine_map/$filename.pip0.5.nonsyn" ]; then
      echo "Downloading and processing $filename"
  
      # Download each file from Google Cloud Storage and save it to data/fine_map
      gsutil cp gs://finngen-public-data-r6/finemapping/summaries/$filename "data/fine_map/$filename"
  
      # Create .pip0.5 filtered files
      # Filtering: remove snps with PIP <= 0.5
      awk 'NR == 1 {print; next} $5 > 0.5' "data/fine_map/$filename" > "data/fine_map/$filename.pip0.5"

      # Create .pip0.5.nonsyn files
      # Fltering: 
      # include only: noncoding fine-mapped loci that did not include any nonsynonymous or splicing variants 
      
      # Currently this means filtering out the following tags: 
      # missense_variant
      # synonymous_variant
      # frameshift_variant
      # splice_acceptor_variant
      # splice_donor_variant
      # stop_gained
      
      # Not sure whether to also filter out: 
      # ? non_coding_transcript_exon_variant
      # ? intergenic_variant
      
      awk '!($14 ~ /missense_variant|synonymous_variant|frameshift_variant|splice_acceptor_variant|splice_donor_variant|stop_gained/)' "data/fine_map/$filename.pip0.5" > "data/fine_map/$filename.pip0.5.nonsyn"

  #awk '$16 != -1' "data/fine_map/finngen_full/$filename" > "data/fine_map/finngen_full/$filename.cs"
  else 
     echo "Skipping downloading + filtering for $filename as already downloaded"
  fi
  
  # if the intermediate files exist delete them
  if [ -f "data/fine_map/$filename" ]; then
  rm "data/fine_map/$filename"
  fi
  
  if [ -f "data/fine_map/$filename.pip0.5" ]; then
  rm "data/fine_map/$filename.pip0.5"
  fi
  
  done < data/fine_map/finngen_snp_filter_list.tsv


################### Step 3: Combine filtered files if they contain at least one snp ##############

echo " Combine the relevant filtered files"
echo " "
echo " "

# Create an empty list to hold filenames
python -c """
from glob import glob
import pandas as pd
pd.concat([
  pd.read_csv(file,sep='\t',header=0) 
  for file in glob('data/fine_map/*nonsyn')
  if pd.read_csv(file,sep='\t',header=0).shape[0] > 0
]).reset_index(drop=True).to_csv('output/fine_map/all_finngen_trait_cs_filtered.tsv', sep='\t', index=False)
"""