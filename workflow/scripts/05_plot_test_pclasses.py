import argparse
import pandas as pd

import stathelper


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
    
    parser.add_argument("result_test_pclasses", type=str,
                        help=("Path to the result file for the test of variants "
                              "by protein classes"))
    parser.add_argument("plot_test_pclasses", type=str,
                        help=("Path to the plot file for the test of variants "
                              "by protein classes"))
    parser.add_argument("plot_title", type=str, help="Title of the plot")
    
    return parser.parse_args(args)




if __name__ == '__main__':
    
    args = parsing()
    
    # Exclude βbridge and π-helix, and sort by the number of protein classes with
    # significant ORs, in the analysis of gnomAD vs ClinVar
    features_test_pclasses = [
        # 'π-helix',
        'Hydrophobic introduced',
        # 'βbridge',
        '310helix',
        # 'Bend',
        'Charge switch',
        'Polar to NonPolar',
        # 'HBondTurn',
        'Big to Small',
        'Aromatic to polar',
        'Small to Big',
        'Aromatic to NonAromatic',
        'Charge lost',
        'Medium-buried (25-50%)',
        'Charge gain',
        'NonPolar to Polar',
        # 'β-sheet/strand',
        'strandβladder',
        'α-helix',
        'Hydrophilic introduced',
        # 'helices',
        'Exposed (>75%)',
        'Core (<5%)',
        'TotalEnergy',
        # 'coils',
        'DisorderpLDDT',
        'OrderpLDDT',
        'Loop',
        'Buried (5-25%)',
        'Contacts',
        'LessTotalEnergy',
        'VanDerWaalsClashes',
        'Medium-exposed (50-75%)',
        'Conserved'
    ]

    
    # All protein classes, sorted by the number of features that have significant
    # ORs, in the analysis of gnomAD vs ClinVar
    pclasses = [
        'C22_transmembrane_signal_receptor',
        'C15_protein_modifying_enzyme',
        'C14.1_hydrolase',
        'C14.6_transferase',
        'C11_gene_specific_transcriptional_regulator',
        'C16_protein_binding_activity_modulator',
        'C23_transporter',
        'C14.5_oxidoreductase',
        'C8_cytoskeletal_protein',
        'C14.3_ligase',
        'C10_extracellular_matrix_protein',
        'C2_RNA_metabolism_protein',
        'C13_membrane_traffic_protein',
        'C4_cell_adhesion_molecule',
        'C17_scaffold_adaptor',
        'C1_DNA_metabolism_protein',
        'C7_chromatin_binding_regulatory_protein',
        'C19_structural_protein',
        'C12_intercellular_signal_molecule',
        'C20_transfer_carrier_protein',
        'C9_defense_immunity_protein',
        'C21_translational_protein',
        'C6_chaperone',
        'C14.4_lyase',
        'C5_cell_junction_protein',
        'C3_calcium_binding_protein',
        'C14.2_isomerase',
        'C24_viral_transposable_element_protein',
        'C18_storage_protein'
    ]
    
    
    results = pd.read_pickle(args.result_test_pclasses)
    
    print(f"Plotting results of variants by protein classes to {args.plot_test_pclasses}...")
    stathelper.plot_heatmap(results, args.plot_title,
                            features_test_pclasses, ['All']+pclasses,
                            figpath=args.plot_test_pclasses)