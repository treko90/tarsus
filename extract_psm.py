#!/usr/bin/env python3
from collections import deque, defaultdict
import csv
import sys

input_file =  sys.argv[1]+'psm.tsv'
output_file = sys.argv[1]+'psm_input4pepquery.tsv'

with open(input_file, 'r', newline='') as fin, open(output_file, 'w', newline='') as fout:
    reader = csv.DictReader(fin, delimiter='\t')
    writer = csv.writer(fout, delimiter='\t')

    for row in reader:
        if float(row['Intensity']) > 0:
            spectrum_str = row['Spectrum']
            d = spectrum_str.split('.')
            file_name = d[0].split('_sub')[0]
            scan_num = str(int(d[1]))
            charge = row['Charge']
            peptide = row['Peptide']
            spectrum_title = f"{file_name}:{scan_num}:{charge}"
            writer.writerow([peptide, spectrum_title])

