import pandas as pd
import os
import argparse

datasets_path = '/sc-scratch/sc-scratch-dh-ukb-intergenics/geno_datasets/'

def count_genes_and_save(name):
    file_path = os.path.join(datasets_path, f'{name}.var.csv')  
    df = pd.read_csv(file_path, sep=',')  
    unique_genes = df['ANN_Gene_ID'].unique() 
    
    output_file = os.path.join(datasets_path+f'{name}/', f'{name}_unique_genes.txt')
    with open(output_file, 'w') as f:
        for gene in unique_genes:
            f.write(f"{gene}\n")
    
    print(f"Unique genes saved to {output_file}")
    return len(unique_genes)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Count and save unique genes in a specified dataset file.")
    parser.add_argument("name", type=str, help="Name of the dataset file (without extension)")
    
    args = parser.parse_args()
    name = args.name
    
    unique_gene_count = count_genes_and_save(name)
    print(f"Number of unique genes in {name}: {unique_gene_count}")
