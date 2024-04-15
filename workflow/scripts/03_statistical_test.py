import argparse
import numpy as np
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
    
    parser.add_argument("path_physicochemical", type=str,
                        help=("Path to the pickle file with all the variants "
                              "and their physicochemical properties"))
    parser.add_argument("path_proteinclass", type=str,
                        help="Path to the protein class pickle file")
    parser.add_argument("result_test_all", type=str,
                        help="Path to the result file for the test of all variants")
    parser.add_argument("plot_test_all", type=str,
                        help="Path to the plot file for the test of all variants")
    parser.add_argument("result_test_pclasses", type=str,
                        help=("Path to the result file for the test of variants "
                              "by protein classes"))
    parser.add_argument("plot_test_pclasses", type=str,
                        help=("Path to the plot file for the test of variants "
                              "by protein classes"))
    
    return parser.parse_args(args)




if __name__ == '__main__':
    
    args = parsing()
    
    # Selected by cecilia in 02_StatisticalTest.ipynb
    features_test_all = [
        'TotalEnergy',
        'LessTotalEnergy',
        'VanDerWaalsClashes',
        'Contacts',
        'Conserved',
        'Core (<5%)',
        'Buried (5-25%)',
        'Medium-buried (25-50%)',
        'Medium-exposed (50-75%)',
        'Exposed (>75%)',
        'DisorderpLDDT',
        'OrderpLDDT',
        'helices',
        'β-sheet/starnd',
        'coils',
        'α-helix',
        'βbridge',
        'strandβladder',
        '310helix',
        'π-helix',
        'HBondTurn',
        'Bend',
        'Loop',
        'Small to Big',
        'Big to Small',
        'Polar to NonPolar',
        'NonPolar to Polar',
        'Hydrophilic introduced',
        'Hydrophobic introduced',
        'Charge switch',
        'Charge lost',
        'Charge gain',
        'Aromatic to NonAromatic',
        'Aromatic to polar',
    ]
    
    # Exclude βbridge and π-helix, and sort by the number of protein classes with
    # significant ORs, in the analysis of gnomAD vs ClinVar
    features_test_pclasses = [
        # 'π-helix',
        'Hydrophobic introduced',
        # 'βbridge',
        '310helix',
        'Bend',
        'Charge switch',
        'Polar to NonPolar',
        'HBondTurn',
        'Big to Small',
        'Aromatic to polar',
        'Small to Big',
        'Aromatic to NonAromatic',
        'Charge lost',
        'Medium-buried (25-50%)',
        'Charge gain',
        'NonPolar to Polar',
        'β-sheet/starnd',
        'strandβladder',
        'α-helix',
        'Hydrophilic introduced',
        'helices',
        'Exposed (>75%)',
        'Core (<5%)',
        'TotalEnergy',
        'coils',
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
    
    
    variants = pd.read_pickle(args.path_physicochemical)
    variants_pclasses = pd.read_pickle(args.path_proteinclass)
    
    # Test of all variants
    sum_results = stathelper.calculate_ODS_mult(features_test_all, variants)
    sum_results = sum_results.dropna(subset=['OR'])
    sum_results = sum_results.sort_values('OR', ascending=False)
    
    print(f"Writing results of all variants to {args.result_test_all}...")
    sum_results.to_pickle(args.result_test_all)
    print(f"Plotting results of all variants to {args.plot_test_all}...")
    stathelper.plot_ODS(sum_results, figpath=args.plot_test_all)
    
    # Test of variants by protein classes
    groupbyName = 'GeneralProteinClass'
    results = stathelper.calculate_ODS_multiple_grouped_by(pclasses, groupbyName,
                                        features_test_pclasses, variants_pclasses)
    
    print(f"Writing results of variants by protein classes to {args.result_test_pclasses}...")
    results.to_pickle(args.result_test_pclasses)
    print(f"Plotting results of variants by protein classes to {args.plot_test_pclasses}...")
    stathelper.plot_heatmap(results, features_test_pclasses, ['All']+pclasses,
                            figpath=args.plot_test_pclasses)