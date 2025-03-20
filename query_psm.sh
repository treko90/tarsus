#!/bin/bash

directory=$1  
fastadir=$2
current_dir=$(pwd)

mkdir -p ${directory}pepquery_logs
mkdir -p ${directory}pepquery_inputs

# copies the pepquery input file and split it into 10-PSM chunks
cp ${directory}second/psm_input4pepquery.tsv ${directory}pepquery_inputs
split ${directory}pepquery_inputs/psm_input4pepquery.tsv -l 10 ${directory}pepquery_inputs/
rm ${directory}pepquery_inputs/psm_input4pepquery.tsv

mv ${fastadir}/protein.faa ${fastadir}/protein.fasta
mv ${fastadir}/protein.fa ${fastadir}/protein.fasta

for file in "${directory}pepquery_inputs"/*; do
	echo "Processing: $file"
	
	filename=$(basename "$file")
	mkdir -p ${directory}pepquery_result/$filename/
	
	# symlinks the fasta file for validation, prevents PepQuery from occasional crashing
	cd $directory
	ln -sf ${current_dir}/${fastadir}/protein.fasta ${current_dir}/${directory}pepquery_result/$filename/sym.fasta
	cd $current_dir
	
	# loads the PepQuery validation with 10-PSM chunk. Modify the command to submit to your computational cluster, add the directory ehre PepQuery is installed
	java -jar /home/labs/pilpel/slavat/pepquery-2.0.2/pepquery-2.0.2.jar -ti 0,1,2,3 -maxLength 50 -tol 20 -itol 0.05 -hc -aa -m 1 -ms ${directory}pepquery_index/ -db ${directory}pepquery_result/$filename/sym.fasta -o ${directory}pepquery_result/$filename/ -e 2 -c 2 -i $file
	
done
