import pandas as pd
import argparse

if __name__ == "__main__":
    file_path_parser = argparse.ArgumentParser()
    file_path_parser.add_argument("--i_path", help="Input path to input dataframe containing variants rsid's and their associated chromosomes.")
    file_path_parser.add_argument("--o_path", help="Output path to save seperated chromosomes files.")
    args = file_path_parser.parse_args()

    file = pd.read_csv(args.i_path, usecols=["ID", "GeneticChromosome"])

    for i in range(1, 23):
        filename = args.o_path + "chr" + str(i) + ".csv"
        
        file[file.GeneticChromosome == i].ID.to_csv(filename, index=False, header=False)

    
# python split_chr.py --i_path "path to csv file" --o_path "output path"