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

# Set working directories and load data
results_folder = '/path/to/results_folder/extended_figs/'
sp_data_folder = '/path/to/spatial_data_folder/C2'

# Create paths and names to results folders for reference regression and cell2location models
ref_run_name = f'{sp_data_folder}/reference_signatures'
run_name = f'{sp_data_folder}/cell2location_map'

t = sc.read_h5ad(run_name + "/sp.h5ad")

# Add 5% quantile, representing confident cell abundance, 'at least this amount is present'
t.obs[t.uns['mod']['factor_names']] = t.obsm['q05_cell_abundance_w_sf']

# Plot figure
plt.clf()
slide = select_slide(t, '5A')
clust_labels = ['C2-7', 'C2-8']
with mpl.rc_context({'figure.figsize': (15, 15)}):
    fig = plot_spatial(
        adata=slide,
        # Labels to show on a plot
        color=clust_labels, labels=clust_labels,
        show_img=True,
        style='fast',
        reorder_cmap=[0, 2],
        # Limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # Size of locations (adjust depending on figure size)
        circle_diameter=6,
        colorbar_position='right'
    )
fig.gca().set_aspect('equal', adjustable='box')
fig.gca().set_xlim(fig.gca().get_xlim()[::-1])
fig.gca().set_ylim(fig.gca().get_ylim()[::-1])
fig.gca().invert_xaxis()

fig
plt.savefig(results_folder + 'humanHYPOMAP_extendedfigure5_oligo.png', dpi=300, bbox_inches='tight')

# Switch to new spatial data folder for C3
sp_data_folder = '/path/to/spatial_data_folder/C3'

ref_run_name = f'{sp_data_folder}/reference_signatures'
run_name = f'{sp_data_folder}/cell2location_map'

t = sc.read_h5ad(run_name + "/sp.h5ad")

# Add 5% quantile, representing confident cell abundance, 'at least this amount is present'
t.obs[t.uns['mod']['factor_names']] = t.obsm['q05_cell_abundance_w_sf']

# Plot figure
plt.clf()
slide = select_slide(t, '2B')
clust_labels = ['C3-13', 'C3-12']
with mpl.rc_context({'figure.figsize': (15, 15)}):
    fig = plot_spatial(
        adata=slide,
        color=clust_labels, labels=clust_labels,
        show_img=True,
        style='fast',
        reorder_cmap=[0, 3],
        # Limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # Size of locations (adjust depending on figure size)
        circle_diameter=6,
        colorbar_position='right'
    )
fig.gca().set_aspect('equal', adjustable='box')
fig.gca().set_xlim(fig.gca().get_xlim()[::-1])
fig.gca().set_ylim(fig.gca().get_ylim()[::-1])
fig.gca().invert_xaxis()

fig
plt.savefig(results_folder + 'humanHYPOMAP_extendedfigure5_ependymal.png', dpi=300, bbox_inches='tight')
