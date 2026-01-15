#!/bin/bash
directory=$1  
fastadir=$2
current_dir=$(pwd)
mkdir -p ${directory}pepquery_logs_base
mv ${fastadir}/protein.faa ${fastadir}/protein.fasta
mv ${fastadir}/protein.fa ${fastadir}/protein.fasta

for file in "${directory}pepquery_input_base"/*; do
	echo "Processing: $file"
	
	filename=$(basename "$file")
	mkdir -p ${directory}pepquery_result_recheck_base/$filename/

	if [ ! -f "${directory}pepquery_result_recheck_base/$filename/psm_rank.txt" ]; then
		cd $directory
		ln -sf ${current_dir}/${fastadir}/protein.fasta ${current_dir}/${directory}pepquery_result_recheck_base/$filename/sym.fasta
		cd $current_dir
		
		# loads the PepQuery validation using the mgf file with the changed precursor masses to match the reference peptide m/z. Modify the command to submit to your computational cluster, add the PepQuery home directory
		java -jar /home/labs/pilpel/slavat/pepquery-2.0.2/pepquery-2.0.2.jar -ti 0,1,2,3 -maxLength 50 -tol 20 -fast -itol 0.01 -fixMod 0 -varMod 1,2 -maxVar 4 -hc -m 1 -s 2 -ms ${directory}mgfall_fakes.mgf -indexType 2 -db ${directory}pepquery_result_recheck_base/$filename/sym.fasta -o ${directory}pepquery_result_recheck_base/$filename/ -e 2 -c 2 -plot -i $file
		fi

done
