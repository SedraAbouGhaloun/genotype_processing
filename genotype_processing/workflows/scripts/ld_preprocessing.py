import pandas as pd
import os
import argparse

def process_dataset(ds_name):
    datasets_path = '/sc-scratch/sc-scratch-dh-ukb-intergenics/geno_datasets/'  # Base path for input datasets
    output_path = f'/sc-scratch/sc-scratch-dh-ukb-intergenics/ld/{ds_name}/'   # Base path for output files

    input_file = f"{datasets_path}{ds_name}/{ds_name}.txt"
    df = pd.read_csv(input_file, sep='\t')

    df = df.rename(columns={
        'chr_name': 'GeneticChromosome',
        'effect_allele': 'REF',
        'other_allele': 'ALT',
        'chr_position': 'POS'
    })

    df = df[['ID', 'POS', 'REF', 'ALT', 'GeneticChromosome']]

    os.makedirs(output_path, exist_ok=True)

    output_file = f'{output_path}{ds_name}.csv'
    print(f"Saving processed file to: {output_file}")
    df.to_csv(output_file, index=False)

    print(f"Processed file saved to: {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process dataset and save reformatted output.")
    parser.add_argument("ds_name", type=str, help="Name of the dataset to process (folder and file prefix)")

    args = parser.parse_args()
    process_dataset(args.ds_name)
