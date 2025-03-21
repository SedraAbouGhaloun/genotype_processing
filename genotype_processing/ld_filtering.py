import os
import pandas as pd
from typing import List 

def process_ld_data(path:str)-> List:
    """
    Loads LD data from multiple ancestries and chromosomes.
    """
    ld = []
    columns = ['CHR_A', 'BP_A', 'SNP_A', 'CHR_B', 'BP_B', 'SNP_B', 'R2']
    ancestries = ['afr', 'amr', 'eur', 'csa', 'eas', 'mid']
    
    for chrom in range(1, 23):
        frames = []
        for ancestry in ancestries:
            file_path = f"{path}/{ancestry}/chr_{chrom}_0.05.ld"
            if os.path.exists(file_path):
                df = pd.read_csv(file_path, sep='\s+')
            else:
                df = pd.DataFrame(columns=columns)
            frames.append(df)
        
        m = pd.concat(frames, ignore_index=True)
        m = m.drop_duplicates(subset=['BP_A', 'SNP_A', 'BP_B', 'SNP_B']).reset_index(drop=True)
        ld.append(m) 
        
    return ld

def get_correlated_variants(ld, variants, max_vars:int=100_000):
    """
    Selects a subset of variants ensuring low LD correlation.
    """
    in_ld_to_something_set = set()
    selected_vars = set()
    
    for var in variants.index:
        if variants['ID'][var] not in in_ld_to_something_set:
            chromosome = variants['GeneticChromosome'][var]
            selected_vars.add(variants['ID'][var])
            in_ld_to_something_set.add(variants['ID'][var])
            ld_set = vars_in_ld_to(variants['ID'][var], chromosome)
            in_ld_to_something_set.update(ld_set)
            
            if len(selected_vars) == max_vars:
                break
    
    return selected_vars


