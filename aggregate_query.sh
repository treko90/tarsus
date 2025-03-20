#!/bin/bash

parent_directory=${1}/  # Change this to your actual directory path

# Define the output file where concatenated data will be stored
output_file="${parent_directory}../psmrankall.txt"

# Clear the output file if it exists, or create a new one
> "$output_file"
# Loop through all directories and subdirectories
find "$parent_directory" -type f -name "psm_rank.txt" | while IFS= read -r file; do
    echo "Processing $file"
    
    # Skip the first line of each psm_rank.txt file and concatenate to output file
    tail -n +2 "$file" >> "$output_file"
    
done

# extract peptide sequences from psmrankall.txt file and save as peptides.txt
awk '{print $1}' $output_file > $parent_directory/../peptides.txt
sed -i '1i\peptide	modification	n	spectrum_title	charge	exp_mass	tol_ppm	tol_da	isotope_error	pep_mass	mz	score	n_db	total_db	n_random	total_random	pvalue	rank	n_ptm	confident	ref_delta_score	mod_delta_score' $output_file


echo "All psm_rank.txt files concatenated into $output_file"

# Define the output file where concatenated data will be stored
output_file="${parent_directory}../ptmall.txt"

# Clear the output file if it exists, or create a new one
> "$output_file"

# Loop through all directories and subdirectories
find "$parent_directory" -type f -name "ptm_detail.txt" | while IFS= read -r file; do
    echo "Processing $file"
    
    # Skip the first line of each psm_rank.txt file and concatenate to output file
    tail -n +2 "$file" >> "$output_file"
    
done

echo "All psm_rank.txt files concatenated into $output_file"

# Define the output file where concatenated data will be stored
output_file="${parent_directory}../mgfall.mgf"

# Clear the output file if it exists, or create a new one
> "$output_file"

# Loop through all directories and subdirectories
find "$parent_directory" -type f -name "psm_rank.mgf" | while IFS= read -r file; do
    echo "Processing $file"
    
    # Skip the first line of each psm_rank.txt file and concatenate to output file
    cat "$file" >> "$output_file"
    
done

echo "All psm_rank.txt files concatenated into $output_file"


mkdir -p ${parent_directory}../pdv_spectra
mkdir -p ${parent_directory}../pdv_spectra_log
		
echo "Processing dataset: $tissue"

# extract fragment ions from the PepQuery outputs
java -jar /home/labs/pilpel/slavat/PDV-2.1.2/PDV-2.1.2.jar -r ${parent_directory}../psmrankall.txt -rt 4 -s ${parent_directory}../mgfall.mgf -st 1 -i ${parent_directory}../peptides.txt -k p -o ${parent_directory}../pdv_spectra/ -ft report
  