#!/bin/bash

fasta_path=$1  
output=$2 
raw_path=$3 
OUTPUT_FILE=${output}first.manifest
dir_path=$(dirname "$fasta_path")
current_dir=$(pwd)

# define the path to your philosopher installation
/home/labs/pilpel/slavat/fragpipe-22.1-build07/tools/Philosopher/philosopher-v5.1.1 workspace --init $dir_path
/home/labs/pilpel/slavat/fragpipe-22.1-build07/tools/Philosopher/philosopher-v5.1.1 database --custom $fasta_path $dir_path
mv *.fas $dir_path
fas_file=$(ls "$dir_path"/*.fas | head -n 1)

# adds new fasta file with decoys to the generic.workflow file. Modify generic.workflow according to your needs
sed "s|database.db-path=/home/labs/pilpel/slavat/domains/fasta/drosophila/GCF_000001215.4/2025-02-18-decoys-protein.faa.fas|database.db-path=${current_dir}/$fas_file|" generic.workflow > ${output}first.workflow

find "$raw_path" -type f -name "*.mzML" | while read -r file; do
	absolute_path=$(realpath "$file")
	echo -e "$absolute_path\t\t\tDDA" >> "$OUTPUT_FILE"
done

# define the path to your fragpipe installation, add cluster submission command
/home/labs/pilpel/slavat/fragpipe-22.1-build07/bin/fragpipe --headless --ram 125 --workflow ${output}first.workflow --manifest $OUTPUT_FILE --workdir $output

mkdir -p ${output}second

python make_peptides.py $output

# define the path to your philosopher installation
/home/labs/pilpel/slavat/fragpipe-22.1-build07/tools/Philosopher/philosopher-v5.1.1 workspace --init $output
/home/labs/pilpel/slavat/fragpipe-22.1-build07/tools/Philosopher/philosopher-v5.1.1 database --custom ${output}pep_candidates.fa 
mv *.fas $output

fas_file=$(ls "$output"/*.fa.fas | head -n 1)

# adds specific peptide search parameters into the second pass workflow file
sed -e "s|database.db-path=|database.db-path=${current_dir}/$fas_file|" -e "s|database.decoy-tag=|database.decoy-tag=rev_|" -e "s|ionquant.minisotopes=2|ionquant.minisotopes=1|" -e "s|ionquant.minscans=3|ionquant.minscans=1|" -e "s|msfragger.clip_nTerm_M=true|msfragger.clip_nTerm_M=false|" -e "s|msfragger.misc.fragger.enzyme-dropdown-1=stricttrypsin|msfragger.misc.fragger.enzyme-dropdown-1=nocleavage|" -e "s|msfragger.remove_precursor_peak=1|msfragger.remove_precursor_peak=0|" -e "s|msfragger.search_enzyme_cut_1=KR|msfragger.search_enzyme_cut_1=@|" -e "s|msfragger.search_enzyme_name_1=stricttrypsin|msfragger.search_enzyme_name_1=nocleavage|" -e "s|msfragger.search_enzyme_nocut_1=|msfragger.search_enzyme_nocut_1=@|" -e "s|msbooster.find-best-rt-model=true|msbooster.find-best-rt-model=false|" ${output}fragpipe-second-pass.workflow > ${output}second.workflow 

# define the path to your fragpipe installation, add cluster submission command
/home/labs/pilpel/slavat/fragpipe-22.1-build07/bin/fragpipe --headless --ram 125 --workflow ${output}second.workflow --manifest ${output}fragpipe-files-second-pass.fp-manifest --workdir ${output}second
