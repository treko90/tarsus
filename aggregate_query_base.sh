#!/bin/bash

parent_directory=${1}/  # Change this to your actual directory path

# Define the output file where concatenated data will be stored
output_file="${parent_directory}../psmrankall_base.txt"

# Clear the output file if it exists, or create a new one
> "$output_file"
# Loop through all directories and subdirectories
find "$parent_directory" -type f -name "psm_rank.txt" | while IFS= read -r file; do
    #echo "Processing $file"
    
    # Skip the first line of each psm_rank.txt file and concatenate to output file
    tail -n +2 "$file" >> "$output_file"
    
done

awk '{print $1}' $output_file > $parent_directory/../peptides.txt
sed -i '1i\peptide	modification	n	spectrum_title	charge	exp_mass	tol_ppm	tol_da	isotope_error	pep_mass	mz	score	n_db	total_db	n_random	total_random	pvalue	rank	n_ptm	confident	ref_delta_score	mod_delta_score' $output_file

echo "All psm_rank.txt files concatenated into $output_file"

# Define the output file where concatenated data will be stored
output_file="${parent_directory}../ptmall_base.txt"

# Clear the output file if it exists, or create a new one
> "$output_file"

# Loop through all directories and subdirectories
find "$parent_directory" -type f -name "ptm_detail.txt" | while IFS= read -r file; do
    #echo "Processing $file"
    
    # Skip the first line of each psm_rank.txt file and concatenate to output file
    tail -n +2 "$file" >> "$output_file"
    
done

echo "All ptm_detail.txt files concatenated into $output_file"

# Define the output file where concatenated data will be stored
output_file="${parent_directory}../mgfall_base.mgf"

# Clear the output file if it exists, or create a new one
> "$output_file"

# Loop through all directories and subdirectories
find "$parent_directory" -type f -name "psm_rank.mgf" | while IFS= read -r file; do
    #echo "Processing $file"
    
    # Skip the first line of each psm_rank.txt file and concatenate to output file
    cat "$file" >> "$output_file"
    
done

echo "All psm_rank.mgf files concatenated into $output_file"
