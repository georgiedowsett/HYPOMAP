###### This script will run the cell2location model using the reference gene signatures from the snRNAseq and the raw ST datasets ######

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
from scipy.sparse import csr_matrix
from cell2location.utils.filtering import filter_genes

# Define paths as variables
results_folder: str = './data/ref_celltype_sigs/240424/C3/'
sp_data_folder: str = './data/input_data/spatial-raw/'
outs_folder: str = './data/c2lmodel-outs/'
barcode_folder: str = './data/input_data/barcodes-to-remove/'
ref_run_name: str = f'{results_folder}/reference_signatures'
run_name: str = f'{results_folder}/cell2location_map'
adata_file: str = f"{ref_run_name}/sc.h5ad"

# Load model and anndata
adata_ref = sc.read_h5ad(adata_file)
mod = cell2location.models.RegressionModel.load(f"{ref_run_name}", adata_ref)

# Export estimated expression in each cluster
def export_estimated_expression(adata_ref: anndata.AnnData) -> pd.DataFrame:
    """
    Export estimated expression in each cluster from the reference anndata object.

    Args:
        adata_ref (anndata.AnnData): The reference anndata object.

    Returns:
        pd.DataFrame: A DataFrame containing the inferred average expression per cluster.
    """
    if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
        inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
    else:
        inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
    inf_aver.columns = adata_ref.uns['mod']['factor_names']
    print(inf_aver.iloc[0:5, 0:5])
    return inf_aver

inf_aver = export_estimated_expression(adata_ref)

def read_and_qc(sample_name: str, barcode: str, path: str = sp_data_folder) -> anndata.AnnData:
    """
    Reads the data for one 10X spatial experiment into the anndata object and calculates QC metrics.
    This function also removes low-quality spots.

    Args:
        sample_name (str): Name of the sample.
        barcode (str): Path to the barcode file.
        path (str): Path to data (default is `sp_data_folder`).

    Returns:
        anndata.AnnData: The processed anndata object.
    """
    barcodes = pd.read_csv(barcode_folder + str(barcode))
    adata = sc.read_visium(path + str(sample_name),
                           count_file='filtered_feature_bc_matrix.h5', load_images=True)
    adata.var_names_make_unique()
    adata.obs['sample'] = sample_name
    adata.var['SYMBOL'] = adata.var_names
    adata.var.rename(columns={'gene_ids': 'ENSEMBL'}, inplace=True)
    adata.var_names = adata.var['SYMBOL']
    adata.var.drop(columns='SYMBOL', inplace=True)
    cond = adata.obs_names.isin(barcodes['Barcode'])
    adata = adata[adata.obs_names[~cond].copy()]

    # Calculate QC metrics
    adata.X = adata.X.toarray()
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    adata.X = csr_matrix(adata.X)
    adata.var['MT'] = [gene.startswith('MT-') for gene in adata.var_names]
    adata.obs['MT_frac'] = adata[:, adata.var['MT'].tolist()].X.sum(1).A.squeeze() / adata.obs['total_counts']

    # Add sample name to obs names
    adata.obs["sample"] = [str(i) for i in adata.obs['sample']]
    adata.obs_names = adata.obs["sample"] + '_' + adata.obs_names
    adata.obs.index.name = 'spot_id'

    return adata

# Read the list of spatial experiments
sample_data = pd.read_csv(sp_data_folder + 'sample_info.csv')

# Read the data into anndata objects
slides = [read_and_qc(i, k, path=sp_data_folder) for i, j, k in zip(sample_data['sample_name'], sample_data['sample_name2'], sample_data['barcodes_to_remove'])]

# Combine anndata objects together
t = slides[0].concatenate(
    slides[1:],
    batch_key="sample",
    uns_merge="unique",
    batch_categories=sample_data['sample_name'],
    index_unique=None
)

# Remove mitochondria-encoded (MT) genes for spatial mapping
t.obsm['MT'] = t[:, t.var['MT'].values].X.toarray()
t = t[:, ~t.var['MT'].values]

# Find shared genes and subset both anndata and reference signatures
intersect = np.intersect1d(t.var_names, inf_aver.index)
t = t[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

# Prepare anndata for cell2location model
cell2location.models.Cell2location.setup_anndata(adata=t, batch_key="sample")

# Create and train the model
mod = cell2location.models.Cell2location(t, cell_state_df=inf_aver, N_cells_per_location=3, detection_alpha=20)
print(mod.view_anndata_setup())

mod.train(max_epochs=30000, batch_size=None, train_size=1, use_gpu=True)

# Plot ELBO loss history during training
plt.clf()
mod.plot_history(1000)
plt.legend(labels=['full data training'])
plt.savefig(results_folder + "model_training")

# Export estimated cell abundance
t = mod.export_posterior(t, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': True})

# Save model and results
mod.save(f"{run_name}", overwrite=True)
t.write(f"{run_name}/sp.h5ad")

# Plot QC
plt.clf()
mod.plot_QC()
plt.savefig(results_folder + "model2_qc")

plt.clf()
fig = mod.plot_spatial_QC_across_batches()
fig.savefig(results_folder + "qc-across-batches")

# Cluster spots into regions using scanpy
sc.pp.neighbors(t, use_rep='q05_cell_abundance_w_sf', n_neighbors=9)
sc.tl.leiden(t, resolution=0.5)
t.obs["region_cluster"] = t.obs["leiden"].astype("category")

# Save anndata object with clusters
t.write(f"{run_name}/sp-clusters.h5ad")

# Export the matrix of values
t.obs[t.uns['mod']['factor_names']] = t.obsm['q05_cell_abundance_w_sf']
x = pd.DataFrame(t.obs)
x.to_csv(results_folder + 'C3_clusters_obs_table.csv')
