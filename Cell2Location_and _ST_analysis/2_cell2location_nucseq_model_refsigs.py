###### This script will create the reference gene expression model from the snRNAseq dataset at each clustering level of the whole human HYPOMAP dataset ######
# Following the cell2location pipeline

# Import packages
import os
import sys
import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import cell2location
from cell2location.models import RegressionModel
import scvi
import matplotlib.pyplot as plt
from cell2location.utils.filtering import filter_genes

# Set environment variables for Theano
os.environ["THEANO_FLAGS"] = 'device=0,floatX=float32,force_device=True'

# Define paths as variables
clusters = ['C1', 'C2', 'C3', 'C4']
base_results_folder: str = "./data/ref_celltype_sigs/240424/"
sp_data_folder: str = './data/input_data/'
nucseq_folder: str = './human_hypo_combined/'

def process_cluster(cluster: str) -> None:
    """
    Processes the given cluster by creating reference gene expression models from snRNAseq data.
    
    Args:
        cluster (str): The cluster identifier (e.g., 'C1', 'C2', etc.)
    """
    print(cluster)

    # Create paths and names to results folders for reference regression and cell2location models
    results_folder = f"{base_results_folder}{cluster}/"
    ref_run_name = f'{results_folder}/reference_signatures'
    run_name = f'{results_folder}/cell2location_map'

    # Load snRNAseq reference data
    adata_ref = anndata.read_h5ad(nucseq_folder + "human_hypo_combined.h5ad")

    # Work on the raw count data
    adata_ref = adata_ref.raw.to_adata()

    # Adding gene names and ENSID
    adata_ref.var = adata_ref.var.rename(columns={"_index": "SYMBOL"})
    adata_ref.var.index = adata_ref.var['SYMBOL']
    del adata_ref.raw

    # Filtering genes
    selected = filter_genes(adata_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.08, nonz_mean_cutoff=1.4)
    adata_ref = adata_ref[:, selected].copy()

    # Prepare anndata for the regression model
    cell2location.models.RegressionModel.setup_anndata(
        adata=adata_ref, 
        batch_key='Sample_ID', 
        labels_key=cluster, 
        categorical_covariate_keys=['Donor_ID', 'sex', 'Dataset']
    )

    mod = RegressionModel(adata_ref)
    print(mod.view_anndata_setup())

    mod.train(max_epochs=250, use_gpu=True)

    # Plot training history
    plt.clf()
    mod.plot_history(20)
    plt.savefig(results_folder + 'elbo-plot-training.png')

    # Export estimated cell abundance
    adata_ref = mod.export_posterior(adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': True})

    # Save model
    mod.save(f"{ref_run_name}", overwrite=True)
    # Save anndata object with results
    adata_file = f"{ref_run_name}/sc.h5ad"
    adata_ref.write(adata_file)

    # Plot model quality control (QC)
    plt.clf()
    mod.plot_QC()
    plt.savefig(results_folder + 'model-qc-plots.png')

def main() -> None:
    """
    Main function to iterate over clusters and process each one.
    """
    for cluster in clusters:
        process_cluster(cluster)

if __name__ == "__main__":
    main()
