import pandas as pd

def exclude_PRS_cohorts(dataframe, target_cohorts, metadata_path):
    """
    Excludes specific cohorts from the dataset based on metadata.
    """
    pgs_metadata = pd.read_csv(metadata_path)
    pgs_metadata['Cohort(s)'] = pgs_metadata['Cohort(s)'].fillna('n')
    
    filtered_pgs_ids = set(pgs_metadata[pgs_metadata['Cohort(s)'].str.contains('|'.join(target_cohorts))]['Polygenic Score (PGS) ID'])
    
    return dataframe[~dataframe['pgs'].str.contains('|'.join(filtered_pgs_ids))].reset_index(drop=True)