import pandas as pd
import argparse
from pathlib import Path
from datetime import datetime

def process_CGD_phenotypes(input_filename, output_filename):
    # provide column names
    column_names = [
            "GENE_SOURCE_ID",
            "Feature Type",
            "Gene Name",
            "CGDID",
            "PMID",
            "Experiment Type",
            "Mutant Type",
            "Allele",
            "Strain background",
            "Phenotype",
            "Chemical",
            "Condition",
            "Details",
            "Reporter",
            "Anatomical Structure",
            "Virulence Model",
            "Species"]
    # Read the .tab file into a DataFrame without header and give column names
    df = pd.read_csv(input_filename, sep='\t', header=None, names=column_names)
    # Split text in the "PMID" column after ":"
    df['PMID'] = df['PMID'].str.split(pat=':', n=1, expand=True)[1]
    # insert a new column called "CGD_REF"
    df.insert(df.columns.get_loc("PMID")+1,"CGD_REF","")
    # Split text in the "PMID" column after "|" and sort into the PMID and CGD_REF columns
    df[['PMID','CGD_REF']] = df['PMID'].str.split(pat='|', n=1, expand=True)
    # Remove "CGD_REF: " from cells in the "CGD_REF" column
    df['CGD_REF'] = df['CGD_REF'].str.replace('CGD_REF: ', '')
    # drop the feature type column
    df.drop(columns=["Feature Type"], inplace=True)
    # Replace empty cells with "no data found"
    df.fillna("no data found", inplace=True)
    # Save the modified DataFrame to a new .tab file
    df.to_csv(output_filename, sep='\t', index=False)
    return

def main():
    parser = argparse.ArgumentParser(
    description="Reformats the CGD phenotype information",
    epilog="written by Helen Davison and Evelina Basenko")
    parser.add_argument('-f','--file', \
                        help="The CGD tab delimited file of phenotype data",
                        required=True)
    args = parser.parse_args()

    # Provide the filenames for the input and output .tab files
    input_filename = args.file
    output_filename = Path(input_filename).stem + "_PROCESSED_"+datetime.today().strftime('%Y-%b-%d') +".tab"
    print("\033[33m{}\033[0;0m".format("\nProcessing file:"))
    print(str(Path(input_filename)))
    process_CGD_phenotypes(input_filename, output_filename)
    print("\033[32m{}\033[0;0m".format('\nNew file saved to:'))
    print(str(Path(output_filename).resolve())+"\n")
    return

main()


