import os
import pandas as pd

def write_dataset_readme(dataset_path, notebook_path, num_genes, contains_all_genes, excluded_cohorts, variant_sources, contains_ld_filter):
    """
    Writes a README file for a dataset, summarizing its details.
    """
    readme_content = f"""
    - Notebook Used for Generation: {notebook_path}
    - Number of Genes: {num_genes}
    - Contains All Genes Condition: {"Yes" if contains_all_genes else "No"}
    - Excluded Cohorts: {", ".join(excluded_cohorts) if excluded_cohorts else "None"}
    - Variant Sources: {variant_sources}
    - Contains LD Filtering: {"Yes" if contains_ld_filter else "No"}
    """
    
    readme_path = os.path.join(os.path.dirname(dataset_path), "README.md")
    with open(readme_path, "w") as readme_file:
        readme_file.write(readme_content.strip() + "\n")

def rename_columns(df):
    return df.rename(columns={"#CHROM": 'chr_name', "POS": "chr_position",
                              "ALT": "effect_allele", "REF": "other_allele",
                              "BETA": "effect_weight"})

def rename_columns_again(df):
    return df.rename(columns={"chr_name": "#CHROM", "chr_position": "POS",
                              "effect_allele": "ALT", "other_allele": "REF",
                              "effect_weight": "BETA"})

__all__ = ["write_dataset_readme"]