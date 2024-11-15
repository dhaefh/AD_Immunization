# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                        AN1792 Project                              -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Written by: Anne Forsyth
# Summary: Estimate cell type abundance with Cell2Location
#
#-------------------------------------------------------------------------------
# Initialization

# Load libraries
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
os.environ["THEANO_FLAGS"] = 'device=cuda,floatX=float32,force_device=True'
import cell2location
from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42
from cell2location.utils.filtering import filter_genes
from cell2location.models import RegressionModel
import torch
from scvi import REGISTRY_KEYS

# Define output folder
output_folder = '/path/to/cell2location/mapping'

# Define path to reference regression output folder
ref_folder = '/path/to/scRNAseq/reference_signatures'

# Load post-regression reference data and Visium data
adata_ref = sc.read(f"{ref_folder}/adata_ref.h5ad")
adata_vis = sc.read('/path/to/all_cohorts_gray.h5ad')

# Define variable for gene names 
adata_vis.var['SYMBOL'] = adata_vis.var.index

# Remove genes prefixed with MT or HB
adata_vis.var['MT_gene'] = [gene.startswith('MT-') for gene in adata_vis.var_names]
adata_vis.var['HB_gene'] = [gene.startswith('HB') for gene in adata_vis.var_names]
adata_vis.obsm['MT'] = adata_vis[:, adata_vis.var['MT_gene'].values].X.toarray()
adata_vis = adata_vis[:, ~adata_vis.var['MT_gene'].values]
adata_vis.obsm['HB'] = adata_vis[:, adata_vis.var['HB_gene'].values].X.toarray()
adata_vis = adata_vis[:, ~adata_vis.var['HB_gene'].values]

# Load reference regression model
mod = cell2location.models.RegressionModel.load(f"{ref_folder}", adata_ref)

# Export reference signatures
if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
  inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}' for i in adata_ref.uns['mod']['factor_names']]].copy()
else:
  inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}' for i in adata_ref.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_ref.uns['mod']['factor_names']
inf_aver.to_csv(f"{ref_folder}/inf_aver.csv")

# Filter reference signatures and spatial data for shared genes
intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()
inf_aver.to_csv(f"{ref_folder}/inf_aver_genes_used.csv")

# Set up cell2location model
cell2location.models.Cell2location.setup_anndata(adata=adata_vis, batch_key="sample_id", continuous_covariate_keys = ["cdr_centered"])
mod = cell2location.models.Cell2location(
  adata_vis, cell_state_df=inf_aver,
  N_cells_per_location=7, 
  detection_alpha=20)
mod.view_anndata_setup()

# Train model in batches 
mod.train(max_epochs=1000, batch_size=2500, train_size=1, accelerator='gpu')

# Plot training curve
mod.plot_history(0)
plt.legend(labels=['full data training'])
plt.savefig(f"{output_folder}/c2l_training_curve.pdf", format = "pdf", bbox_inches = "tight")
plt.close()
plt.cla()
plt.clf()

# Export summary of the posterior distribution of cell type abundance (computing quantiles directly)
adata_vis = mod.export_posterior(
  adata_vis, use_quantiles=True,
  add_to_obsm=["q05", "q50", "q95", "q0001"],
  sample_kwargs={'batch_size': 2500, 'use_gpu': True})

# Save model
mod.save(f"{output_folder}", overwrite=True)

# Save spatial anndata object
adata_file = f"{output_folder}/adata_vis.h5ad"
adata_vis.write(adata_file)

# Generate reconstruction plot
mod.plot_QC(summary_name = "q50")
plt.savefig(f"{output_folder}/reconstruction_plot.pdf", format = "pdf", bbox_inches = "tight")
plt.close()
plt.cla()
plt.clf()

