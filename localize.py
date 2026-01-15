import os
import polars as pl
import sys
from pyteomics import mgf
import numpy as np
from polars import col
import re
from pyascore import PyAscore

def localize_shift(unbased_df,unfiltered_df,depquery_spectra_path,prec):
    totest = unbased_df.rename({'peptide':'dependent'}).join(unfiltered_df[['dependent','base','origin','mistranslatedPos','destination']].unique(),on='dependent')
    #totest = unfiltered_df
    #spectra = mgf.IndexedMGF(os.getcwd()+'/revisions/hoMS/GAPDH/2d/mzml/dependent_spectra.mgf') # gapdh
    #spectra = mgf.IndexedMGF(os.getcwd()+'/revisions/hoMS/hemoglobin/2d/fragpipe/mgfall_dep.mgf') #hb 
    spectra = mgf.IndexedMGF(depquery_spectra_path)
    
    aa_mono = {
        "A": 71.03711,
        "R": 156.10111,
        "N": 114.04293,
        "D": 115.02694,
        "C": 103.00919,
        "E": 129.04259,
        "Q": 128.05858,
        "G": 57.02146,
        "H": 137.05891,
        "I": 113.08406,
        "L": 113.08406,
        "K": 128.09496,
        "M": 131.04049,
        "F": 147.06841,
        "P": 97.05276,
        "S": 87.03203,
        "T": 101.04768,
        "W": 186.07931,
        "Y": 163.06333,
        "V": 99.06841,
    }
    
    substitution_shifts = {
        aa_from: {
            aa_to: round(aa_mono[aa_to] - aa_mono[aa_from], 4)
            for aa_to in aa_mono if aa_to != aa_from
        }
        for aa_from in aa_mono
    }
    
    def get_var_mod_positions(
        best_seq: str,
        peptide: str,
        mod_mass: float,
        aux_mod_pos: np.ndarray,
        aux_mod_mass: np.ndarray,
        mass_tol: float = prec,
    ):
        """
        Return 1-based residue positions that carry *your* variable mod
        (defined by mod_mass), accounting for aux_mod_pos/aux_mod_mass.
        """
        # 1) Build a map of static (aux) mass per position
        static_mass = {}
        for pos, mass in zip(aux_mod_pos, aux_mod_mass):
            pos = int(pos)
            if pos < 0:
                continue
            # 0 is N-term; handle separately if you allow N-term variable mods
            static_mass[pos] = static_mass.get(pos, 0.0) + float(mass)
    
        positions = []
        aa_pos = 0
        i = 0
    
        while i < len(best_seq):
            c = best_seq[i]
    
            if c.isalpha():
                aa_pos += 1
                # check if this residue is followed by a [mass]
                if i + 1 < len(best_seq) and best_seq[i + 1] == "[":
                    j = best_seq.find("]", i + 2)
                    if j == -1:
                        raise ValueError("Unbalanced brackets in best_sequence")
                    mass_str = best_seq[i + 2 : j] if best_seq[i+2] in "+-" else best_seq[i+1 : j]
                    total_mass = float(mass_str)
    
                    # subtract known static (aux) mass at this position
                    static = static_mass.get(aa_pos, 0.0)
                    var_mass = total_mass - static
    
                    # if the residual mass matches your mod_mass -> this is your PTM
                    if abs(var_mass - mod_mass) <= mass_tol:
                        positions.append(aa_pos)
    
                    i = j + 1
                else:
                    i += 1
            else:
                i += 1
    
        return positions
    
    WATER_MASS  = 18.01056         # H2O added for neutral peptide
    PROTON_MASS = 1.007276466812  # mass of a proton
    
    
    def mz_diff(peptide: str, observed_mz: float, charge: int):
        """
        Compute theoretical monoisotopic m/z for an unmodified peptide
        and the difference vs an observed m/z.
    
        Args:
            peptide: peptide sequence (one-letter AA codes, no mods)
            observed_mz: observed precursor m/z
            charge: precursor charge state (e.g. 2, 3, 4)
    
        Returns:
            (theoretical_mz, delta_da, delta_ppm)
            where:
              theoretical_mz = calculated m/z
              delta_da       = observed_mz - theoretical_mz (in Da)
              delta_ppm      = (delta_da / theoretical_mz) * 1e6
        """
        # neutral monoisotopic mass = sum(residues) + H2O
        try:
            neutral_mass = WATER_MASS + sum(AA_MONO[aa] for aa in peptide)
        except KeyError as e:
            raise ValueError(f"Unknown amino acid in sequence: {e.args[0]}")
    
        theoretical_mz = (neutral_mass + charge * PROTON_MASS) / charge
        delta_mz = observed_mz - theoretical_mz
        delta_ppm = (delta_mz / theoretical_mz) * 1e6
        delta_neutral = delta_mz * charge   # <-- this is the ~71 Da
    
        return theoretical_mz, delta_mz, delta_ppm, delta_neutral
    
    def ascore_to_probability(ascore):
        P_wrong = 10 ** (-ascore / 10)
        P_correct = 1 - P_wrong
        return P_wrong
    
    sure_pos,ascore_list,ascore_dict = [],[],dict()
    for row in totest.iter_rows(named=True):
    
        positions = [int(pos) for pos in re.findall(r'@(\d+)', row['modification'])]
        masses = [float(mass) for mass in re.findall(r'\[([\d.]+)\]', row['modification'])]
        mod_group="ACDEFGHIKLMNPQRSTVWY"
        spec = spectra.get_spectrum(row['spectrum_title'])
        #mod_mass = mz_diff(row['base'],row['observed_mz'],row['charge'])[3]
        mod_mass = substitution_shifts[row['origin']][row['destination']]
        peptide = row['base']
        aux_mod_pos=np.array(positions, dtype=np.uint32)
        aux_mod_mass=np.array(masses, dtype=np.float32)
        
        mz    = spec["m/z array"]
        inten = spec["intensity array"]
        z     = spec["params"].get("charge", [2])[0]  # e.g. [2]
        
        ascore = PyAscore(
            bin_size=100.0,
            n_top=10,
            mod_group=mod_group,
            mod_mass=mod_mass,
     
            mz_error=prec,
            fragment_types="by"
        )
        
        ascore.score(
            mz_arr=mz,
            int_arr=inten,
            peptide=peptide,
            n_of_mod=1,
            max_fragment_charge=z-1,
            #max_fragment_charge=1,
            aux_mod_pos=aux_mod_pos,
            aux_mod_mass=aux_mod_mass,
        )
    
        # 1) Map signature indices -> peptide positions of modifiable residues
        modifiable_positions = [
            i for i, aa in enumerate(peptide, start=1) if aa in mod_group
        ]
         #len(modifiable_positions) == len(signature)
        
        # 2) Take best localization
        best = ascore.pep_scores[0]
        
        
        sig  = best["signature"]
        
        # All modifiable sites with the PTM (for n_of_mod > 1)
        modified_modifiable_indices = np.where(sig == 1)[0]

        #print(ascore.ascores)
        #print(ascore.pep_scores[0])
        #print(ascore.alt_sites)
        #print(row)
        #print()
        # Map back to 1-based peptide positions
        best_peptide_positions = [modifiable_positions[i] for i in modified_modifiable_indices]
        if len(best_peptide_positions) != 1: print(best_peptide_positions)
        if len(best_peptide_positions) == 1 and best_peptide_positions[0] == row['mistranslatedPos']+1 and ascore.pep_scores[0]['weighted_score'] > ascore.pep_scores[1]['weighted_score']:
        #print(best_peptide_positions[0]==row['mistranslatedPos']+1,len(best_peptide_positions))
    
            #print(ascore.best_sequence,row['dependent'],row['mistranslatedPos'],mod_mass)
            #sure_pos.append(alltissue_dict[row['spectrum_title']])
            sure_pos.append(row['spectrum_title'])
    
            ascore_list.append(ascore_to_probability(ascore.ascores[0]))
            ascore_dict[row['spectrum_title']] = ascore.ascores[0]
            
        #if row['spectrum_title'] == '01294_H01_P013197_S00_N08_R1_sub:21952:2':
            #rint(ascore.pep_scores)
            #print(alltissue_dict[row['spectrum_title']])
    
            #print(ascore.best_score)
            #print(substitution_shifts[row['origin']][row['destination']])
            #print(ascore.ascores)
            #print(ascore.pep_scores)
            #print(ascore.alt_sites)
            #print(row)
            #print(ascore.best_sequence,row['dependent'],row['mistranslatedPos'],mod_mass)
        #print(len(sure_pos))
    return unfiltered_df.filter(col('Spectrum').is_in(sure_pos))


ipc = pl.read_ipc(sys.argv[1]+'subs.ipc')

deps = pl.read_csv(sys.argv[1]+'psmrankall.txt',separator='\t')
base = pl.read_csv(sys.argv[1]+'psmrankall_base.txt', separator='\t').with_columns(
    pl.col("spectrum_title").str.replace(r"_[^_]*$", "")
)

nonbasehits = deps.join(base.group_by(['peptide','spectrum_title']).agg(pl.col('score').max()),on=['spectrum_title']).with_columns((pl.col('score')-pl.col('score_right')).alias('delta_score')).filter(pl.col('delta_score') > 0)
ipc_localized = localize_shift(nonbasehits,ipc,sys.argv[1]+'mgfall.mgf',0.01).with_columns(pl.col('Protein_right').str.split('.').list.first()).rename({'Protein_right':'ensembl'})

print(ipc_localized)
ipc_localized.write_ipc(sys.argv[1]+'subs_localized.ipc')