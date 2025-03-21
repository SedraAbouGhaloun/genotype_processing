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


def vars_in_ld_to(var, chromosome, ld):
    chr_ld = ld[chromosome - 1]
    chr_ld_df = chr_ld[(chr_ld['SNP_A'] == var) | (chr_ld['SNP_B'] == var)]
    LDs_set = set(chr_ld_df.SNP_B)
    LDs_set.update(set(chr_ld_df.SNP_A))

    return LDs_set


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



def get_remaining_ld(ld, df1, df2, r2_threshold=0.05):
    """
    Extracts LD pairs not fully within df1 or df2.

    Parameters:
    - ld (pd.DataFrame): LD DataFrame with columns SNP_A, SNP_B, R2.
    - df1 (pd.DataFrame): First dataset with column 'ID'.
    - df2 (pd.DataFrame): Second dataset with column 'ID'.
    - r2_threshold (float): Minimum R² threshold for filtering LD pairs.

    Returns:
    - pd.DataFrame: Remaining LD pairs not in ld_within_df1 or ld_within_df2.
    """
    df1_ids = set(df1['ID'])
    df2_ids = set(df2['ID'])
    
    ld = ld[ld['R2'] > r2_threshold]
    ld = ld[(ld['SNP_A'].isin(df1_ids | df2_ids)) & (ld['SNP_B'].isin(df1_ids | df2_ids))]

    ld_within_df1 = ld[
        (ld['SNP_A'].isin(df1_ids)) & 
        (ld['SNP_B'].isin(df1_ids)) & 
        (~ld['SNP_A'].isin(df2_ids)) & 
        (~ld['SNP_B'].isin(df2_ids))
    ]

    ld_within_df2 = ld[
        (ld['SNP_A'].isin(df2_ids)) & 
        (ld['SNP_B'].isin(df2_ids)) & 
        (~ld['SNP_A'].isin(df1_ids)) & 
        (~ld['SNP_B'].isin(df1_ids))
    ]

    ld_rest = ld[~ld.index.isin(ld_within_df1.index) & ~ld.index.isin(ld_within_df2.index)]

    return ld_rest

def modified_jaccard_index(df1, df2, ld_matrix, LD_threshold=0.05):
    variants_df1 = set(df1['ID'])
    variants_df2 = set(df2['ID'])
    shared_variants = variants_df1.intersection(variants_df2)
    len_shared_variants = len(shared_variants)
    # print(f'number of shared variants: {str(len_shared_variants)}')
    
    ld_variants_df1 = set(ld_matrix[ld_matrix['R2'] >= LD_threshold]['SNP_A'])
    ld_variants_df2 = set(ld_matrix[ld_matrix['R2'] >= LD_threshold]['SNP_B'])
    
    ld_in_union = ld_variants_df1.union(ld_variants_df2)
    len_ld_in_union = len(ld_in_union)
    
    intersection_ld_shared = shared_variants.intersection(ld_in_union)
    len_intersection_ld_shared = len(intersection_ld_shared)
    # print(f'number of shared varints in LD: {str(len_intersection_ld_shared)}')

    
    total_variants_union = variants_df1.union(variants_df2)
    jaccard_index = (len_shared_variants + len_ld_in_union - len_intersection_ld_shared) / len(total_variants_union)
    
    return jaccard_index



