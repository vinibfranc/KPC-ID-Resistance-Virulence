import sys
import os
import glob
import csv
import pandas as pd

REPORT_FOLDER = "results/kraken2/report/"
PATHOGENS_LIST = "lista_patogenos.csv"
os.system("mkdir -p results/kraken_genus_species_strains/")
os.system("mkdir -p results/kraken_pathogens/")

def filter_kraken_report():
    
    patho_file = open(PATHOGENS_LIST, 'r')
    reader = csv.reader(patho_file)
    patho_dict = {}
    for row in reader:
        taxid, specie = row
        patho_dict[taxid] = specie
    # print(patho_dict)

    patho_df = pd.DataFrame(patho_dict.items())
    # print(patho_df)

    patho_file.close()

    files = glob.glob(REPORT_FOLDER+'*.kreport')
    for file in files:
        # print(file)
        kreport = pd.read_csv(file, sep='\t', lineterminator='\n')
        kreport.iloc[:, 5] = kreport.iloc[:, 5].str.strip()

        not_species = kreport.loc[(kreport.iloc[:, 3] != "S") & (kreport.iloc[:, 3] != "S1") & (kreport.iloc[:, 3] != "G")].index
        kreport.drop(not_species, inplace=True)

        format_seq = os.path.basename(file)
        #print(format_seq)

        print("Species reported by Kraken2 on ", format_seq)
        print(kreport)

        kreport.to_csv(r'results/kraken_genus_species_strains/'+format_seq, sep='\t', index=None, header=True)

        # Compare 1th column of dayaframe (list of pathogens) and 4th column of dataframe (species reported), excluding all lines where there are none match (no CNS pathogen, just microorganism)
        print("Pathogens identified: ")

        kreport.iloc[:, 4] = kreport.iloc[:, 4].astype(int)
        patho_df.iloc[:, 0] = patho_df.iloc[:, 0].astype(int)

        pathogens_identified = pd.merge(kreport, patho_df, how='inner', left_on=kreport.iloc[:, 4], right_on=patho_df.iloc[:, 0])
        pathogens_identified.drop(pathogens_identified.columns[[0, 7, 8]], axis=1, inplace=True)
        print(pathogens_identified)

        pathogens_identified.to_csv(r'results/kraken_pathogens/'+format_seq, sep='\t', index=None, header=True)
    
def main():
    filter_kraken_report()

if __name__ == '__main__':
    main()