# filter gwas data to genomic regions shown in figure


##################### Figure 4a. #################

# 6q15
# 93100001-99500000 hg19 ucsc browser
# 92390283-99052124 liftover from hg19 to hg38
# in hg38 88000001-93100000 ucsc browser

output_dir="output/gwas_ss_filt"

fig_4a_data=("$output_dir/ra_uk_bb.h.filt.tsv"
             "$output_dir/t1d_uk_bb.h.filt.tsv"
              "$output_dir/hypo_uk_bb.h.filt.tsv"
              "$output_dir/finngen_atopic_derm.filt.tsv"
             )

echo "Filtering Figure 4a datasets"
for i in "${!fig_4a_data[@]}"; do
    gwas_sum_stat="${fig_4a_data[$i]}"

    echo " "
    echo "$gwas_sum_stat"

    if [[ "$gwas_sum_stat" == *"uk"* ]]; then

        awk -F'\t' 'NR == 1 || $3 == 6 && $4 >= 88000001 && $4 <= 93100000' "$gwas_sum_stat" >  "${gwas_sum_stat%.tsv}.fig4a.tsv"

    elif [[ "$gwas_sum_stat" == *"finngen"* ]]; then

        awk -F'\t' 'NR == 1 || $1 == 6 && $2 >= 88000001 && $2 <= 93100000' "$gwas_sum_stat" >  "${gwas_sum_stat%.tsv}.fig4a.tsv"

    fi
done




######################## Figure 4 b. #######################
# 4p15.2

# 21300001-27700000 hg19
# 21298378-27698378 hg38
echo " "
echo " "
echo "Filtering Figure 4b datasets"

fig_4b_data=("$output_dir/ra_uk_bb.h.filt.tsv"
             "$output_dir/t1d_uk_bb.h.filt.tsv"
             )

for i in "${!fig_4b_data[@]}"; do
    gwas_sum_stat="${fig_4b_data[$i]}"

    echo " "
    echo "$gwas_sum_stat"

    awk -F'\t' 'NR == 1 || $3 == 4 && $4 >= 21298378 && $4 <= 27698378' "$gwas_sum_stat" >  "${gwas_sum_stat%.tsv}.fig4b.tsv"

done

##################### Figure 4c ######################### 
# 11p14.1
# hg38 liftover 27178454-30978453
# hg19 27200001-31000000
fig_4c_data=("$output_dir/endo_uk_bb.h.filt.tsv"
             "$output_dir/ovary_cys_uk_bb.h.filt.tsv"
             "$output_dir/menorrhagia_uk_bb.h.filt.tsv"
             "$output_dir/age_meno_uk_bb.h.filt.tsv"
             )

echo " "
echo " "
echo "Filtering figure 4 c datasets"

for i in "${!fig_4c_data[@]}"; do
    gwas_sum_stat="${fig_4c_data[$i]}"
    
    echo " "
    echo "$gwas_sum_stat"
    
    if [[ "$gwas_sum_stat" == *"age_meno"* ]]; then #menopause files have a different file structure
        
awk -F: 'NR > 1 {
    # Split the first column into chrom, pos, ref, alt
    split($1, fields, ":");
    # Print the new columns (chrom, pos, ref, alt) followed by the remaining columns starting from $2
    # We will ensure that there are no extra spaces between columns by using printf correctly
    printf "%s\t%s\t%s\t%s", fields[1], fields[2], fields[3], fields[4];
    # Print the rest of the columns starting from $2
    for (i=2; i<=NF; i++) {
        printf "\t%s", $i;  # Ensures that each field is separated by a single tab
    }
    print "";  # Ensure a newline is printed at the end of each row
}
' "$gwas_sum_stat" > "${gwas_sum_stat%.tsv}.fig4c.noheader.tsv"

# Now do figure 4 filtering

awk -F'\t' '$1 == 11 && $5 >= 27178454 && $5 <= 30978453'  "${gwas_sum_stat%.tsv}.fig4c.noheader.tsv" >  "${gwas_sum_stat%.tsv}.fig4c.tsv"

# After processing rows, add the header manually
{
    echo -e "chrom\tb1\tb2\tb3\tpos\tref\talt\tminor_allele\tminor_AF\tlow_confidence_variant\tn_complete_samples\tAC\tytx\tbeta\tse\ttstat\tpval";
    cat "${gwas_sum_stat%.tsv}.fig4c.noheader.tsv";
} > "${gwas_sum_stat%.tsv}.fig4c.tsv"

# # Clean up the temporary file# rm "${gwas_sum_stat%.tsv}.fig4c.noheader.tsv"

    else 
    
       awk -F'\t' 'NR == 1 || $3 == 11 && $4 >= 27178454 && $4 <= 30978453' "$gwas_sum_stat" >  "${gwas_sum_stat%.tsv}.fig4c.tsv"
    fi
    
done 