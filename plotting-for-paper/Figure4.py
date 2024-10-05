# Load packages 
import os
import sys
import scanpy as sc
import anndata
import pandas as pd
import numpy as np
os.environ["THEANO_FLAGS"] = 'device=cuda,floatX=float32,force_device=True'
import cell2location
from cell2location.models import RegressionModel
import scvi
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import rcParams
from cell2location.utils import select_slide
from scipy import sparse
from cell2location.utils.filtering import filter_genes
from cell2location.plt import plot_spatial

# Set directories
results_folder = '/path/to/results_folder/figures/'
sp_data_folder = '/path/to/spatial_data_folder/'

# Create paths and names to results folders for reference regression and cell2location models
ref_run_name = f'{sp_data_folder}/reference_signatures'
run_name = f'{sp_data_folder}/cell2location_map'

# Load data 
t = sc.read_h5ad(run_name + "/sp.h5ad")
# Add 5% quantile, representing confident cell abundance, 'at least this amount is present',
t.obs[t.uns['mod']['factor_names']] = t.obsm['q05_cell_abundance_w_sf']

# Plot figure 
plt.clf()
slide = select_slide(t, '2B')
clust_labels = ['C4-373', 'C4-374', 'C4-375']
with mpl.rc_context({'figure.figsize': (15, 15)}):
    fig = plot_spatial(adata=slide, color=clust_labels, labels=clust_labels, show_img=True, style='fast', max_color_quantile=0.992, circle_diameter=6, colorbar_position='right')
plt.savefig(results_folder + 'humanHYPOMAP_POMC_spatial.png', dpi=300, bbox_inches='tight')
