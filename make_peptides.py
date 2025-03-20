import polars as pl
import pandas as pd
import sys
import os

print(os.getcwd()+'/{}'.format(sys.argv[1]))

def single_character_substitutions(original: str, target: str, substitute: str):
	"""
	Generate all possible single character substitutions in the string and 
	return the substitutions along with the positions of each substitution.
	param original: The original string (e.g., "banana").
	:param target: The character to be replaced (e.g., 'n').
	:param substitute: The character to substitute with (e.g., 's').
	:return: A tuple of two lists: (list of substitutions, list of substitution positions).
	"""
	substitutions = []
	positions = []
	# Find all indices where the target character appears
	indices = [i for i, char in enumerate(original) if char == target]
	   
	# Replace the target character at each index with the substitute
	for index in indices:
		# Create a new string with the substitution
		new_string = original[:index] + substitute + original[index + 1:]
		substitutions.append(new_string)
		positions.append(index) 
	return substitutions, positions


subs_allowed = pd.read_csv('subs_exclude.csv',index_col=0).rename(index={'I/L': 'I'}).to_dict()
peptides = pl.read_csv(os.getcwd()+'/{}peptide.tsv'.format(sys.argv[1]),separator='\t')
peptides = peptides.filter(pl.col('Intensity') > 0)

with open(os.getcwd()+'/{}/pep_candidates.fa'.format(sys.argv[1]),mode='w') as candidates:
	#print(peptides['Peptide'])
	peptides = set(peptides['Peptide'])
	toprint = []
	#print(peptides)
	for i in peptides:
		for origin in subs_allowed:
			for dest in subs_allowed[origin]:
				if subs_allowed[origin][dest] == False and origin != dest:
					if origin not in i:
						#toprint.append('>'+i+'_n'+'\n'+i)
						pass
					else:
						variants,positions = single_character_substitutions(i,origin,dest)
						#toprint.append('>'+i+'_b'+'\n'+i)
						#print(variants,positions)
						for var,pos in zip(variants,positions):
							if var not in peptides:
								toprint.append('>'+i+'_{}{}{}'.format(origin,pos,dest)+'\n'+var)
	for pep in toprint:
		print(pep,file=candidates)

