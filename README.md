# Workflow to calculate feature enrichment on pathogenic vs benign protein variants

With this Snakemake workflow you can compute the Fisher exact test to calculate feature enrichment on multiple datasets. You can provide one or more benign and one or more pathogenic datasets, and the workflow will compute the Fisher exact test for every combination of benign and pathogenic datasets.

## Installation and usage

1. Clone the repository.

2. Install and activate the conda environment specified on `environment.yaml`.

3. Specify the datasets to analyze in the `config/config.yml` file.

4. Run the snakemake workflow:

```bash
# Run the workflow in dry mode to check that everything is set up correctly
snakemake -np

# Submite the workflow to the cluster
sbatch out_ibex/run_snakemake.sh
```

## Output

The workflow will generate several directories in `results/`, and their contents are as follows:

- `results/03_statistical_test`: Pickle files with the results of the Fisher exact test for every combination of benign and pathogenic datasets.
- `results/04_plot_all`: Plots of the Odds Ratios with their respective CIs for every analysis.
- `results/05_plot_pclasses`: Heatmaps of the Odds Ratios for the tests of the variants grouped by protein class.