
# module load CBI
# module load miniforge3/24.7.1-0
# module load R/4.4.2
# conda activate bmi-206-group
# #conda install conda-forge::gsutil


mkdir -p output/fine_map

########### Part 1: Determine: what finngen files to dowload ##################

# a. Get all the file names for the finemapping data v6

# get a full list of avaliable snp-filter trait files
# snp-filter susie datasets are 
# snps that are part of good quality credible sets as reported in PHENONAME.cred.summary.tsv file.
gsutil ls gs://finngen-public-data-r6/finemapping/summaries/ |\
           grep '\.SUSIE\.snp\.filter\.tsv$' |\
           sed 's|.*/||' > data/fine_map/finngen_snp_filter_list.tsv


# b. the traits tested/file names don't exactly match what's avaliable in the bucket for the manifest file- 
# so let's overlap them
gsutil cp gs://finngen-public-data-r6/summary_stats/R6_manifest.tsv "data/fine_map/R6_manifest.tsv"
        
  
# we overlap them while doing filtering step 1: filter for traits with cases > 1,000         
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

matching_files = data.frame(matching_files)

data.table::fwrite(matching_files,
                    "data/fine_map/finngen_snp_filter_list.tsv",
                    col.names=F
                  )
'

########### Step 2. Download the files, and do basic processing ##############

#module load htslib/1.21

# download + basic process all the relevant files
while read -r filename; do

  if [ ! -f "data/fine_map/$filename.pip0.5.nonsyn" ]; then
  
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
  fi
  
rm "data/fine_map/$filename"
rm "data/fine_map/$filename.pip0.5"
  
  done < data/fine_map/finngen_snp_filter_list.tsv


################### Step 3: Combine filtered files if they contain at least one snp ##############


# Create an empty list to hold filenames
output_file="output/fine_map/all_finngen_trait_cs_filtered.tsv"

> "$output_file"  # Initialize the output file (clear it if it exists)

# Loop through each .snp.filter.tsv file in the folder

for file in data/fine_map/*.pip0.5.nonsyn; do
  # Check if the file has more than 1 line
  line_count=$(wc -l < "$file")
  
  if [ "$line_count" -gt 1 ]; then
      cat "$file" >> "$output_file"
      echo "" >> "$output_file"  # Add a newline after each file for separation
  fi
done
