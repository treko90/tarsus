#!/bin/bash

# Define the directory peptide.txt 
results_dir=$2  # Replace with your actual directory path
files_dir=$1

mkdir -p ${results_dir}pepquery_index

/home/labs/pilpel/slavat/pepquery-2.0.2/pepquery-2.0.2.jar index -r -i $files_dir -o ${results_dir}pepquery_index
