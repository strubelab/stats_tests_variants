import argparse
import numpy as np
import pandas as pd

DESCRIPTION = ("")

def parsing(args: list=None) -> argparse.Namespace:
    """
    Creates the argument parser instance and applies it to the command line
    input

    Args:
        args (list, optional): List of the arguments to be parsed (only to be
            used for testing). If none is provided, it is taken from sys.argv.
            Defaults to None.

    Returns:
        argparse.Namespace
    """

    parser = argparse.ArgumentParser(description=DESCRIPTION)
    
    parser.add_argument("benign_dataset", type=str,
                        help="Path to the benign dataset")
    parser.add_argument("pathogenic_dataset", type=str,
                        help="Path to the pathogenic dataset")
    parser.add_argument("destination", type=str,
                        help="Path to the destination pickle file")
    
    return parser.parse_args(args)


if __name__ == '__main__':

    args = parsing()
    
    benign = pd.read_pickle(args.benign_dataset)
    benign['Source'] = 'Benign'

    pathogenic = pd.read_pickle(args.pathogenic_dataset)
    pathogenic['Source'] = 'Pathogenic'

    concatenated_df = pd.concat([benign, pathogenic], ignore_index=True)
    
    # Remove duplicated rows that are both in benign and pathogenic, keeping
    # only one the benign ones
    # This is assuming that there are no duplicates within the datasets, only
    # between them
    select = (concatenated_df.duplicated(subset=['#CHROM', 'POS', 'REF', 'ALT'], keep=False) 
              & (concatenated_df['Source'] == 'Benign'))
    index_to_remove = concatenated_df[select].index
    concatenated_df = concatenated_df.drop(index_to_remove).reset_index(drop=True)
    
    # Split the amino acids into two columns
    concatenated_df[['REFAA', 'ALTAA']] = (concatenated_df['Amino_acids']
                                           .str.split('/', expand=True))

    # Replace concatenated_df.secondary_structure == 'nan' with 'L'
    concatenated_df.secondary_structure = (concatenated_df.secondary_structure
                                           .replace(np.nan, 'L'))
    
    # Normalize residue accessibility
    # Dictionary mapping amino acids to their respective x values
    x_values = {
        'Ala': 129.0,
        'Arg': 274.0,
        'Asn': 195.0,
        'Asp': 193.0,
        'Cys': 167.0,
        'Glu': 223.0,
        'Gln': 225.0,
        'Gly': 104.0,
        'His': 224.0,
        'Ile': 197.0,
        'Leu': 201.0,
        'Lys': 236.0,
        'Met': 224.0,
        'Phe': 240.0,
        'Pro': 159.0,
        'Ser': 155.0,
        'Thr': 172.0,
        'Trp': 285.0,
        'Tyr': 263.0,
        'Val': 174.0
    }

    # Create the 'ACCESIBILITY_NORMALIZED' column
    concatenated_df['ACCESIBILITY_NORMALIZED'] = concatenated_df.apply(
            lambda row: row['accessibility'] / x_values[row['Residue']], axis=1)
    
    # Write to pickle
    print(f"Writing to {args.destination}...")
    concatenated_df.to_pickle(args.destination)
