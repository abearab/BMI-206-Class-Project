#!/bin/bash

#$ -S /bin/bash     # run job as a Bash shell [IMPORTANT]
#$ -cwd             # run job in the current working directory

# Use this email address:
#$ -M isobel.beasley@ucsf.edu

# Send yourself an email when the job:
# aborts abnormally (fails), begins, and ends successfully
#$ -m abe

# Multithreaded (SMP) job: must run on one node 
#$ -pe smp 1

# The name of the job:
#$ -N gwas_down

# How much RAM per slot
#$ -l mem_free=10G

#$ -l scratch=2G      # job requires up to 2 GiB of local /scratch space

# The maximum running time of the job in hours:mins:sec (converted to 10 minutes):
#$ -l h_rt=0:20:00

module load CBI 
#module load miniforge3/24.7.1-0
#conda activate bmi-206-group
module load htslib/1.21

# Download relevant gwas summary statistics


# Define arrays for URLs, output file names, and disease names
urls=(
     "https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/M05.gwas.imputed_v3.both_sexes.tsv.bgz"
#    "ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90077001-GCST90078000/GCST90077873/harmonised/34662886-GCST90077873-EFO_0000685.h.tsv.gz"

    "ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90014001-GCST90015000/GCST90014023/harmonised/34012112-GCST90014023-EFO_0001359.h.tsv.gz"
    
    "https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/HYPOTHYROIDISM.gwas.imputed_v3.both_sexes.tsv.bgz"
    #"ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90077001-GCST90078000/GCST90077731/harmonised/34662886-GCST90077731-EFO_1001055.h.tsv.gz"
    
    "https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/N80.gwas.imputed_v3.both_sexes.tsv.bgz"
    "ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90077001-GCST90078000/GCST90077821/harmonised/34662886-GCST90077821-EFO_0001065.h.tsv.gz"
    
    "https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/20002_1349.gwas.imputed_v3.female.tsv.bgz"
    "ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90077001-GCST90078000/GCST90077793/harmonised/34662886-GCST90077793-HP_0000138.h.tsv.gz"
    
    "https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/20002_1556.gwas.imputed_v3.female.tsv.bgz"
    #"ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90077001-GCST90078000/GCST90077933/harmonised/34662886-GCST90077933-HP_0000132.h.tsv.gz"
   
    "https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/3581_irnt.gwas.imputed_v3.both_sexes.tsv.bgz"
     #"ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90079001-GCST90080000/GCST90079064/harmonised/GCST90079064.h.tsv.gz"
     
    "https://storage.googleapis.com/finngen-public-data-r6/summary_stats/finngen_R6_ATOPIC_STRICT.gz"
)

data_dir="data/gwas_ss"
output_dir="output/gwas_ss_filt"

# Add corresponding output files here
output_files=(
              "$data_dir/ra_uk_bb.h.tsv"
              "$data_dir/t1d_uk_bb.h.tsv"
              "$data_dir/hypo_uk_bb.h.tsv"
              "$data_dir/endo_uk_bb.h.tsv"
              "$data_dir/ovary_cys_uk_bb.h.tsv"
              "$data_dir/menorrhagia_uk_bb.h.tsv"
              "$data_dir/age_meno_uk_bb.h.tsv"
              "$data_dir/finngen_atopic_derm"
              )


 # Add corresponding disease names here                  
disease_names=(
              "Rheumatoid Arthritis"
              "Type 1 Diabetes"
              "Hypothyroidism"
              "Endometriosis"
              "Ovary Cysts"
              "Menorrhagia"
              "Age of Menopause"
              "Atopic Dermatitis"
              ) 

# Loop over each element in the arrays
for i in "${!urls[@]}"; do
    url="${urls[$i]}"
    output_file="${output_files[$i]}"
    disease_name="${disease_names[$i]}"

    # Download and unzip the data
    echo " "
    echo "Now processing summary statistics for ${disease_name}"
    echo "output_file: ${output_file}"
    echo "url: ${output_file}"
    echo " "
    echo " "

    if [ ! -f "${output_file}.gz" ]; then
        wget -nc -c -q -O "${output_file}.gz" "$url"
    
       if [[ "$url" == *"broad"* &&  ! -f "${output_file}.bgz" ]]; then
           echo "Downloading data for ${disease_name}..."
           wget -nc -c -q -O "${output_file}.bgz" "$url"
       fi
    fi
    
    if [ ! -f "$output_file" ]; then
       echo " "
       echo " "
       echo "Unzip data for ${disease_name}..."
       if [[ "$url" == *"broad"* ]]; then
       bgzip -d "${output_file}.bgz"
       else
       gunzip "${output_file}.gz"
       fi
    fi
    
    filt_file="$(echo "$output_file" | sed "s|^$data_dir|$output_dir|" | sed 's|\.tsv$|.filt.tsv|')"
    
    if [[ "$output_file" == *"finngen"* ]]; then
       filt_file="$filt_file.filt.tsv"
    fi
    
    echo " "
    echo "Downloaded file saved as $output_file"
    echo "Filtered file saved as $filt_file"
    echo " "
    
    if [[ "$output_file" == *"uk"* ]]; then
        if [[ "$output_file" == *"broad"* ]]; then #files from broad have different structure
            awk ' NR == 1 || $1 ~ /^6/ || $1 ~ /^4/  || $1 ~ /^11/ ' "$output_file" > "$filt_file"

            #awk -F'\t' 'NR == 1 || $1 == 6 || $1 == 4 || $1 == 11' "$output_file" > "$filt_file"
        else 
        # Filtering and counting based on chromosome and position for "ukbb" files
            awk -F'\t' 'NR == 1 || $3 == 6 || $3 == 4 || $3 == 11' "$output_file" > "$filt_file"
            
        fi
        
        # Print statistics
        echo "${disease_name} number of sites tested for chromosome 4, 6, 11 "
        "$filt_file" | wc -l

        # echo "${disease_name} number of sites for chromosome 4"
        # awk -F'\t' 'NR == 1 || $3 == 4' "$filt_file" | wc -l
        # 
        # echo "${disease_name} number of sites for chromosome 11"
        # awk -F'\t' 'NR == 1 || $3 == 11' "$filt_file" | wc -l
        
    elif [[ "$output_file" == *"finngen"* ]]; then
    
        # Filtering and counting based on chromosome and position for "finngen" files
        awk -F'\t' 'NR == 1 || $1 == 6 || $1 == 4 || $1 == 11' "$output_file" > "$filt_file"
        
        # Print statistics
        echo "${disease_name} number of sites tested for chromosome 4, 6, 11 "
        "$filt_file" | wc -l
        # echo "${disease_name} number of sites tested for chromosome 6"
        # awk -F'\t' 'NR == 1 || $1 == 6' "$filt_file" | wc -l
        # 
        # echo "${disease_name} number of sites for chromosome 4"
        # awk -F'\t' 'NR == 1 || $1 == 4' "$filt_file" | wc -l
        # 
        # echo "${disease_name} number of sites for chromosome 11"
        # awk -F'\t' 'NR == 1 || $1 == 11' "$filt_file" | wc -l
    fi  

    echo "Processing complete for ${disease_name}."
done

# Original scrappy code ... 


# # Rheumatoid Arthritis
# # GCST90077873 - 331754 European
# wget -O data/ra_uk_bb.h.tsv.gz ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90077001-GCST90078000/GCST90077873/harmonised/34662886-GCST90077873-EFO_0000685.h.tsv.gz
# gunzip data/ra_uk_bb.h.tsv.gz
# # awk -F'\t' '$3 == 6 || $3 == 4' ra_uk_bb.h.tsv > ra_uk_bb_filt.h.tsv
# # awk -F'\t' '$3 == 6 && $4 >= 91673441 && $4 <= 92230175' ra_uk_bb.h.tsv | wc -l
# echo "RA number of sites tested for chromosome 6"
# awk -F'\t' '$3 == 6' ra_uk_bb.h.tsv | wc -l
# echo "RA number of sites for chromosome 4"
# awk -F'\t' '$3 == 4' ra_uk_bb.h.tsv | wc -l
# awk -F'\t' '$3 == 6 || $3 == 4 || $3 == 11' data/ra_uk_bb.h.tsv > data/ra_uk_bb_filt.h.tsv
# 
# 
# # Type 1 diabetes 
# # n = 520580 
# wget -O data/t1d_uk_bb.h.tsv.gz ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90014001-GCST90015000/GCST90014023/harmonised/34012112-GCST90014023-EFO_0001359.h.tsv.gz
# gunzip t1d_uk_bb.h.tsv.gz
# 
# echo  "T1D number of sites tested for chromosome 5"
# awk -F'\t' '$3 == 6' t1d_uk_bb.h.tsv.gz | wc -l
# awk -F'\t' '$3 == 4' t1d_uk_bb.h.tsv.gz | wc -l
# 
# # Hypothyroidism 
# # n= 329052
# wget -O hypo_uk_bb.h.tsv.gz ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90077001-GCST90078000/GCST90077731/harmonised/34662886-GCST90077731-EFO_1001055.h.tsv.gz
# gunzip hypo_uk_bb.h.tsv.gz 
# awk -F'\t' '$3 == 6' hypo_uk_bb.h.tsv  | wc -l
# 
# # Endometro 331754
# wget -O endo_uk_bb.h.tsv.gz ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90077001-GCST90078000/GCST90077821/harmonised/34662886-GCST90077821-EFO_0001065.h.tsv.gz
# gunzip endo_uk_bb.h.tsv.gz
# awk -F'\t' '$3 == 11' endo_uk_bb.h.tsv  | wc -l
# awk -F'\t' '$3 == 11' menorrhagia_uk_bb.h.tsv  | wc -l
# awk -F'\t' '$3 == 11' ovaryc_uk_bb.h.tsv  | wc -l
# awk -F'\t' '$1 == 11' age_meno_uk_bb.h.tsv  | wc -l
# # menorrhagia  331754
# wget -O menorrhagia_uk_bb.h.tsv.gz ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90077001-GCST90078000/GCST90077933/harmonised/34662886-GCST90077933-HP_0000132.h.tsv.gz
# gunzip menorrhagia_uk_bb.h.tsv.gz
# 
# wget -O https://storage.googleapis.com/finngen-public-data-r6/summary_stats/finngen_R6_M13_RHEUMA.gz 
# 
# #ovarian cyst 331754
# wget -O ovaryc_uk_bb.h.tsv.gz ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90077001-GCST90078000/GCST90077793/harmonised/34662886-GCST90077793-HP_0000138.h.tsv.gz
# gunzip ovaryc_uk_bb.h.tsv.gz
# 
# # age of menopause 140172
# wget -O age_meno_uk_bb.h.tsv.gz ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90079001-GCST90080000/GCST90079064/harmonised/GCST90079064.h.tsv.gz
# gunzip age_meno_uk_bb.h.tsv.gz
