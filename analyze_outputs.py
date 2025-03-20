import polars as pl
import os
import re
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import sys

aa_mapping = {"Ala": "A","Arg": "R","Asn": "N","Asp": "D","Cys": "C","Gln": "Q","Glu": "E","Gly": "G","His": "H","Xle": "I","Leu": "I","Lys": "K","Met": "M","Phe": "F","Pro": "P","Ser": "S","Thr": "T","Trp": "W","Tyr": "Y","Val": "V"}


def extract_pdv_ions(input_string):
    """
    Extract values starting with 'b' or 'y' followed by numbers, removing any '+' signs
    and returning unique values in a list.
    
    Args:
        input_string (str): Input string containing values like 'b3', 'y4+++', etc.
        
    Returns:
        list: List of unique values with just letter and number
    """
    # Remove brackets and split by comma
    if type(input_string) == str:
        items = input_string.strip('[]').split(', ')
        
        # Use list comprehension to extract letter and number
        values = [re.match(r'([by])(\d+)', item).group(0) 
                  for item in items 
                  if re.match(r'([by])(\d+)', item)]
        
        # Remove duplicates using set and convert back to list
        return sorted(list(set(values)))
    else:
        return -1

# parses the ptmall.txt and extracts PSM not overlapping with either PTM or alternative peptide sequences
def save_singles(sample):
    ptm = pl.read_csv(os.getcwd()+'/{}/ptmall.txt'.format(sys.argv[1]),separator='\t',new_columns=['peptide','modification','n','spectrum_title','charge','exp_mass','tol_ppm','tol_da','isotope_error','pep_mass','mz','score','n_db','total_db','n_random','total_random','pvalue','rank','ptm_spectrum_title','ptm_peptide','ptm_charge','ptm_exp_mass','ptm_pep_mass','ptm_tol_ppm','ptm_tol_da','ptm_isotope_error','ptm_modification','ptm_score'],truncate_ragged_lines=True,ignore_errors=True)
    
    # puts origin and destination and position columns
    pattern = r'([A-Za-z]{3})->([A-Za-z]{3}).*@(\d+)\['
    ptm = ptm.with_columns([
    pl.col("ptm_modification").str.extract(pattern, 1).alias("origin"),
    pl.col("ptm_modification").str.extract(pattern, 2).alias("destination"),
    pl.col("ptm_modification").str.extract(pattern, 3).cast(pl.Int64).alias("mistranslatedPos")])
    
    ptm = ptm.with_columns([(pl.col("mistranslatedPos") - 1).alias("mistranslatedPos")])
    ptm = ptm.with_columns([pl.col("origin").replace(aa_mapping).alias("origin"),pl.col("destination").replace(aa_mapping).alias("destination")])
    
    # select only such PSMs which were assigned to SAVs
    ptm_subonly = ptm.filter(pl.col("spectrum_title").is_in(ptm.group_by("spectrum_title").agg(all_contains = pl.col("ptm_modification").str.contains("->").all()).filter(pl.col("all_contains"))["spectrum_title"]))
    
    # select only valid range of mistranslated positions, removes weird AAs
    df_valid = ptm_subonly.filter(
    (pl.col("mistranslatedPos") >= 0) &
    (pl.col("mistranslatedPos") < pl.col("ptm_peptide").str.len_bytes()))
    
    # creates a modified peptide sequence according to pepquery output
    df_valid = df_valid.with_columns(
    (
        pl.col("ptm_peptide").str.slice(0, pl.col("mistranslatedPos"))
        +
        pl.col("destination")
        +
        pl.col("ptm_peptide").str.slice(pl.col("mistranslatedPos") + 1)
    ).alias("modified_peptide"))
    
    # corrects for I/L ambiguity
    df = df_valid.with_columns([
    pl.col("peptide").str.replace_all(r"[IL]", "J").alias("peptide_mod"),
    pl.col("modified_peptide").str.replace_all(r"[IL]", "J").alias("modified_peptide_mod")])

    group_stats = df.group_by("spectrum_title").agg([pl.col("ptm_score").max().alias("max_ptm_score"),pl.col("ptm_score").filter(pl.col("ptm_score") == pl.col("ptm_score").max()).count().alias("max_count")])
    second_max = df.group_by("spectrum_title").agg([pl.col("ptm_score").filter(pl.col("ptm_score") < pl.col("ptm_score").max()).max().alias("second_max_ptm_score")])
    df_merged = df.join(group_stats, on="spectrum_title", how="left").join(second_max, on="spectrum_title", how="left").with_columns([pl.col('second_max_ptm_score').fill_null(1)])

    df_selected = df_merged.filter((pl.col("max_count") == 1) & (pl.col("ptm_score") == pl.col("max_ptm_score")))

    sure = df_selected.filter((pl.col("peptide_mod") == pl.col("modified_peptide_mod")))

    
    # take only PSMs with one suggested SAV
    sure = sure.group_by('spectrum_title').len().filter((pl.col('len') == 1)).join(sure,on='spectrum_title')

    return sure


def extract_base(df):
    # extracts the reference peptide sequence
    names = []
    for i in df.iter_rows(named=True):
        names.append(i['Protein'].split('_')[0])
    return df.with_columns(pl.Series("base", names))

def generate_ion_series(df):
    """
    Generate b and y ion series in a single column.
    Y ions range from len(dependent)-position to len(dependent)
    """
    return df.with_columns([
        pl.concat_list([
            pl.Series(name='ion_series', values=[
                [f'b{i}' for i in range(pos+1, len(dep) + 1)] + 
                [f'y{i}' for i in range(len(dep)-pos+2, len(dep)+1)]
                for pos, dep in zip(df['position'], df['dependent'])
            ])
        ]).alias('ion_series')
    ])

def mistranslations(sample):
    pdvexport = pl.read_csv(os.getcwd()+'/{}/pdv_spectra/PDVExportResult.txt'.format(sys.argv[1]),separator='\t')
    psmdict = {i:dict() for i in set(pdvexport['peptide'])}
    
    # extract ionseries from the tested PSMs
    for psm in pdvexport.iter_rows(named=True):
        ionseries = extract_pdv_ions(psm['ions'])
        if ionseries != -1:
            psmdict[psm['peptide']][psm['Spectrum_Title']] = ionseries
            
    first = pl.read_csv(os.getcwd()+'/{}/psmrankall.txt'.format(sys.argv[1]),separator='\t')
    first = first.filter(pl.col('confident') == 'Yes') # extract PSM's confidently validated by PepQuery
    firstset = set(first['spectrum_title'])
    
    # look into ptmall.txt and "save" the peptides unambigously assigned to the alternative sequence
    savedset = set(save_singles(sample)['spectrum_title'])
    
    # extract reference peptide (first search) PSM's, maximum intensity PSM represents the quantity
    first_psms = pl.read_csv(os.getcwd()+'/{}/psm.tsv'.format(sys.argv[1]),separator='\t')
    first_psms = first_psms.rename({'Intensity':'first Intensity'})
    first_psms = first_psms.rename({'Peptide': "base"}).filter(pl.col('first Intensity') > 0).with_columns(pl.lit(sys.argv[1].split('/')[-1]).alias("sample"))
    first_psms = first_psms.with_columns([pl.max('first Intensity').over('base','Charge').alias('Intensity_max')]).filter(pl.col('first Intensity') == pl.col('Intensity_max')).drop('Intensity_max')
    
    # extract the alternative peptides (second search) PSM's, maximum intensity PSM represents the quantity
    second_all = extract_base(pl.read_csv(os.getcwd()+'/{}/second/psm.tsv'.format(sys.argv[1]),separator='\t'))
    second_all = second_all.rename({'Peptide':'dependent','Intensity':'second Intensity'})
    second_all = second_all.filter(pl.col('second Intensity') > 0).with_columns(pl.lit(sys.argv[1].split('/')[-1]).alias("sample"))
    second_all = second_all.with_columns()
    
    new = []
    for i in second_all['Spectrum']:
        i = i.split('.')
        new.append(i[0].split('_sub')[0]+':'+i[1].lstrip('0')+':'+i[-1])
    second_all = second_all.with_columns(pl.Series("Spectrum", new))
    second_psms = second_all.filter((pl.col('Spectrum').is_in(firstset)) | (pl.col('Spectrum').is_in(savedset)))
    
    # for each peptide create a list of expected ions
    expected_ions = generate_ion_series(second_psms[['dependent','Protein']].unique().with_columns([
            pl.col('Protein')
            .str.extract(r'\D(\d+)\D*$', group_index=1)  # More flexible pattern
            .cast(pl.Int64)
            .alias('position')
        ]))[['dependent','ion_series']].to_pandas().set_index('dependent').to_dict()['ion_series']
    
    # retain only PSM's with two fragments covering the substituted position
    twofragment_psms = []
    for peptide in expected_ions:
        if peptide in psmdict:
            for psm in psmdict[peptide]:
                if len(set(psmdict[peptide][psm]).intersection(set((expected_ions[peptide])))) > 1:
                    twofragment_psms.append(psm)
    second_psms = second_psms.filter(pl.col('Spectrum').is_in(twofragment_psms))
    
    # only equally chrged reference and alternative peptide precursors are paired 
    second_psms = second_psms.with_columns([pl.max('second Intensity').over('dependent','Charge').alias('Intensity_max')]).filter(pl.col('second Intensity') == pl.col('Intensity_max')).drop('Intensity_max')
    
    # connect reference and their derived alternative peptides, 
    # filter out the N-terminal substitutions and assign origin, destination amino acids 
    # and sub position in the protein, calculate log10 abundance ratios and save the substitution table as ipc in the results directory
    joinedsearches = second_psms.join(first_psms[['Spectrum','base','Prev AA','Next AA','Protein','first Intensity','Retention','Charge','Protein Start','Protein End','Mapped Proteins','Hyperscore','Extended Peptide','Number of Missed Cleavages','Is Unique']],on=['base','Charge'])
    joinedsearches = joinedsearches.with_columns(pl.col("Protein").str.extract(r"_(?:\D*)(\d+)", 1).cast(pl.Int64).alias("mistranslatedPos")).filter(pl.col('mistranslatedPos') > 0)
    joinedsearches = joinedsearches.with_columns([
            pl.col("Protein").str.extract(r'_(\w)\d+(\w)', 1).alias("origin"),
            pl.col("Protein").str.extract(r'_(\w)\d+(\w)', 2).alias("destination"),
            (pl.col("Protein Start_right")+pl.col('mistranslatedPos')-1).alias("misProt")])
    joinedsearches = joinedsearches.with_columns(np.log10(pl.col('second Intensity')/pl.col('first Intensity')).alias('dpbp'))
    joinedsearches.write_ipc(os.getcwd()+'/'+sys.argv[1]+'subs.ipc')
    return (joinedsearches,first_psms)


print(mistranslations(sys.argv[1]))