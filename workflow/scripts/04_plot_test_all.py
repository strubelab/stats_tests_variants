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
    
    parser.add_argument("result_test_all", type=str,
                        help="Path to the file with the test results of all variants")
    parser.add_argument("plot_test_all", type=str,
                        help="Path to the plot file for the test of all variants")
    parser.add_argument("plot_title", type=str, help="Title of the plot")
    
    return parser.parse_args(args)




if __name__ == '__main__':
    
    args = parsing()
    
    sum_results = pd.read_pickle(args.result_test_all)
    
    print(f"Plotting results of all variants to {args.plot_test_all}...")
    stathelper.plot_ODS(sum_results, title=args.plot_title,
                        figpath=args.plot_test_all)
