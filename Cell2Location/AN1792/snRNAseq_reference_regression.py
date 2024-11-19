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
# Summary: Estimate snRNAseq cell type signatures with Cell2Location
#
#-------------------------------------------------------------------------------

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
ref_folder = '/path/to/snRNAseq/reference_signatures'

# Load snRNAseq reference anndata object
ref = sc.read_h5ad('/path/to/snRNAseq/ref.h5ad')

# Define variable for gene names 
ref.var['SYMBOL'] = ref.var.index

# Define gene filtering parameters for gene filter plot
adata = ref.copy()
cell_count_cutoff=5
cell_percentage_cutoff2=0.03
nonz_mean_cutoff=1.12

# Generate gene filter plot using source code
# https://cell2location.readthedocs.io/en/latest/_modules/cell2location/utils/filtering.html#filter_genes
adata.var["n_cells"] = np.array((adata.X > 0).sum(0)).flatten()
adata.var["nonz_mean"] = np.array(adata.X.sum(0)).flatten() / adata.var["n_cells"]

cell_count_cutoff = np.log10(cell_count_cutoff)
cell_count_cutoff2 = np.log10(adata.shape[0] * cell_percentage_cutoff2)
nonz_mean_cutoff = np.log10(nonz_mean_cutoff)

gene_selection = (np.array(np.log10(adata.var["n_cells"]) > cell_count_cutoff2)) | (
  np.array(np.log10(adata.var["n_cells"]) > cell_count_cutoff)
  & np.array(np.log10(adata.var["nonz_mean"]) > nonz_mean_cutoff)
)
gene_selection = adata.var_names[gene_selection]
adata_shape = adata[:, gene_selection].shape

fig, ax = plt.subplots()
ax.hist2d(
  np.log10(adata.var["nonz_mean"]),
  np.log10(adata.var["n_cells"]),
  bins=100,
  norm=mpl.colors.LogNorm(),
  range=[[0, 0.5], [1, 4.5]],
)
ax.axvspan(0, nonz_mean_cutoff, ymin=0.0, ymax=(cell_count_cutoff2 - 1) / 3.5, color="darkorange", alpha=0.3)
ax.axvspan(
  nonz_mean_cutoff,
  np.max(np.log10(adata.var["nonz_mean"])),
  ymin=0.0,
  ymax=(cell_count_cutoff - 1) / 3.5,
  color="darkorange",
  alpha=0.3,
)

plt.vlines(nonz_mean_cutoff, cell_count_cutoff, cell_count_cutoff2, color="darkorange")
plt.hlines(cell_count_cutoff, nonz_mean_cutoff, 1, color="darkorange")
plt.hlines(cell_count_cutoff2, 0, nonz_mean_cutoff, color="darkorange")
plt.xlabel("Mean non-zero expression level of gene (log)")
plt.ylabel("Number of cells expressing gene (log)")
plt.title(f"Gene filter: {adata_shape[0]} cells x {adata_shape[1]} genes")

# Save plot
fig.savefig(f"{ref_folder}/filter_genes.pdf", format = "pdf", bbox_inches = "tight")
plt.close()
plt.cla()
plt.clf()

# Filter genes in reference object
selected = filter_genes(ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)
print(f"{len(selected)} genes selected")
ref = ref[:, selected].copy()

# Set up regression model 
cell2location.models.RegressionModel.setup_anndata(adata=ref,
                                                   batch_key='batch',
                                                   labels_key='ct')
mod = RegressionModel(ref)
mod.view_anndata_setup() 

# Train model  
mod.train(max_epochs=500, accelerator='gpu', batch_size = 2500, train_size = 1, lr=0.002) 

# Export summary of the posterior distribution of expression signatures
ref = mod.export_posterior(
  ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': True}
)

# Save model
mod.save(ref_folder, overwrite=True)

# Save reference anndata object
adata_file = f"{ref_folder}/adata_ref.h5ad"
ref.write(adata_file)

# Generate reconstruction plots using source code (replaces mod.plot_QC - issue saving plots when running on Quest)
# https://cell2location.readthedocs.io/en/latest/_modules/cell2location/models/reference/_reference_model.html#RegressionModel.plot_QC

# Define default reconstruction plot parameters
summary_name = "means"
use_n_obs = 1000
scale_average_detection = True

# First plot
super(type(mod), mod).plot_QC(summary_name = summary_name, use_n_obs = use_n_obs)
plt.savefig(f"{ref_folder}/first_reconstruction_plot.pdf", format = "pdf", bbox_inches = "tight")
plt.close()
plt.cla()
plt.clf()

# Second plot
inf_aver = mod.samples[f"post_sample_{summary_name}"]["per_cluster_mu_fg"].T
if scale_average_detection and ("detection_y_c" in list(mod.samples[f"post_sample_{summary_name}"].keys())):
  inf_aver = inf_aver * mod.samples[f"post_sample_{summary_name}"]["detection_y_c"].mean()
aver = mod._compute_cluster_averages(key=REGISTRY_KEYS.LABELS_KEY)
aver = aver[mod.factor_names_]

plt.hist2d(
  np.log10(aver.values.flatten() + 1),
  np.log10(inf_aver.flatten() + 1),
  bins=50,
  norm=mpl.colors.LogNorm())

plt.xlabel("Mean expression for every gene in every cluster")
plt.ylabel("Estimated expression for every gene in every cluster")

plt.savefig(f"{ref_folder}/second_reconstruction_plot.pdf", format = "pdf", bbox_inches = "tight")
plt.close()
plt.cla()
plt.clf()

# Plot training curve
mod.plot_history(0)
plt.savefig(f"{ref_folder}/reference_training_curve.pdf", format = "pdf", bbox_inches = "tight")
plt.close()
plt.cla()
plt.clf()
