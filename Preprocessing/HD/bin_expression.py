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
# Summary: Map Space Ranger output to segmented nuclei and generate plots to evaluate segmentation
#
#-------------------------------------------------------------------------------

# Load libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import anndata
import geopandas as gpd
import scanpy as sc
from tifffile import imread, imwrite
from csbdeep.utils import normalize
from stardist.models import StarDist2D
from shapely.geometry import Polygon, Point
from scipy import sparse
from matplotlib.colors import ListedColormap
import os
from stardist import render_label
import pickle
import math

# Define general output folder
output_folder = "/path/to/preprocessing/output/folder/"

# Define paths to DAPI images + Space Ranger output
img_folder = '/path/to/DAPI/images/'
spaceranger_folder = '/path/to/spaceranger/output/'

# Define list of samples
samples = ['B9', 'A9', 'A34995E']

# Bin expression for each sample
for sample in samples:
  
  # Define sample subfolder
  result_folder = f'{output_folder}{sample}'
    
  # Load polys
  file = open(f"{result_folder}/polys.pkl","rb")
  data = file.read()
  polys = pickle.loads(data)
  file.close()
    
  # Initialize list to store Polygon geometries
  geometries = []

  # Iterate through each nucleus
  for nuclei in range(len(polys['coord'])):
    
    # Extract coordinates for the current nuclei and convert them to (y, x) format
    coords = [(y, x) for x, y in zip(polys['coord'][nuclei][0], polys['coord'][nuclei][1])]
    
    # Create a Polygon geometry 
    geometries.append(Polygon(coords))

  # Create a GeoDataFrame using the Polygon geometries
  gdf = gpd.GeoDataFrame(geometry=geometries)
    
  # Append sample to nuclei IDs
  gdf['id'] = [f"{sample}_{i+1}" for i, _ in enumerate(gdf.index)]
    
  # Load Space Ranger output and initialize anndata object
  in_names = [sample in name for name in os.listdir(spaceranger_folder)]
  index = [in_names.index(i) for i in in_names if i == True] 
  name = os.listdir(spaceranger_folder)[index[0]]
  raw_h5_file = f'{spaceranger_folder}{name}/outs/binned_outputs/square_002um/filtered_feature_bc_matrix.h5'
  adata = sc.read_10x_h5(raw_h5_file)
    
  # Load spatial coordinates
  tissue_position_file = f'{spaceranger_folder}{name}/outs/binned_outputs/square_002um/spatial/tissue_positions.parquet'
  df_tissue_positions=pd.read_parquet(tissue_position_file)
  df_tissue_positions = df_tissue_positions.set_index('barcode')
  df_tissue_positions['index']=df_tissue_positions.index
    
  # Add spatial coordinates to anndata object meta data
  adata.obs =  pd.merge(adata.obs, df_tissue_positions, left_index=True, right_index=True)
    
  # Create a GeoDataFrame from the coordinates
  geometry = [Point(xy) for xy in zip(df_tissue_positions['pxl_col_in_fullres'], df_tissue_positions['pxl_row_in_fullres'])]
  gdf_coordinates = gpd.GeoDataFrame(df_tissue_positions, geometry=geometry)
    
  # Check which coordinates are in a cell nucleus - this returns multiple rows per barcode if barcode appears in multiple nuclei
  # NA values are assigned to barcodes not in a cell nucleus
  result_spatial_join = gpd.sjoin(gdf_coordinates, gdf, how='left', predicate='within')
    
  # Identify nuclei associated barcodes and barcodes that are not in overlapping nuclei
  result_spatial_join['is_within_polygon'] = ~result_spatial_join['index_right'].isna()
  barcodes_in_overlapping_polygons = pd.unique(result_spatial_join[result_spatial_join.duplicated(subset=['index'])]['index'])
  result_spatial_join['is_not_in_an_polygon_overlap'] = ~result_spatial_join['index'].isin(barcodes_in_overlapping_polygons)
    
  # Filter for barcodes in a cell nucleus and not in overlapping nuclei
  barcodes_in_one_polygon = result_spatial_join[result_spatial_join['is_within_polygon'] & result_spatial_join['is_not_in_an_polygon_overlap']]
  filtered_obs_mask = adata.obs_names.isin(barcodes_in_one_polygon['index'])
  filtered_adata = adata[filtered_obs_mask,:]

  # Add the results of the point spatial join to the anndata object
  filtered_adata.obs =  pd.merge(filtered_adata.obs, barcodes_in_one_polygon[['index','geometry','id','is_within_polygon','is_not_in_an_polygon_overlap']], left_index=True, right_index=True)
    
  # Group the data by unique nucleus IDs
  groupby_object = filtered_adata.obs.groupby(['id'], observed=True)
    
  # Extract the raw expression from the anndata object
  counts = filtered_adata.X

  # Identify the number of unique nuclei and total number of genes
  N_groups = groupby_object.ngroups
  N_genes = counts.shape[1]
    
  # Initialize a sparse matrix to store the summed gene counts for each nucleus
  summed_counts = sparse.lil_matrix((N_groups, N_genes))

  # Iterate over each unique polygon to calculate the sum of gene counts
  polygon_id = []
  row = 0
  for polygons, idx_ in groupby_object.indices.items():
    summed_counts[row] = counts[idx_].sum(0)
    row += 1
    polygon_id.append(polygons)
    
  # Create and anndata object from the summed count matrix
  summed_counts = summed_counts.tocsr()
  grouped_filtered_adata = anndata.AnnData(X=summed_counts,obs=pd.DataFrame(polygon_id,columns=['id'],index=polygon_id),var=filtered_adata.var)
    
  # Store the area of each nucleus in the GeoDataframe
  gdf['area'] = gdf['geometry'].area
    
  # Delete geometry variable from filtered_adata
  del filtered_adata.obs['geometry']
    
  # Save anndata objects
  filtered_adata.write(f"{result_folder}/filtered_adata.h5ad")
  grouped_filtered_adata.write(f"{result_folder}/binned_adata.h5ad")
    
  # Save gdf with area info for all polygons
  if os.path.exists(f"{result_folder}/gdf") == False:
    os.mkdir(f"{result_folder}/gdf")
  gdf.to_file(f"{result_folder}/gdf/gdf.shp")
    
  # Save areas to csv
  areas = pd.DataFrame(gdf[['id','area']])
  areas.to_csv(f"{result_folder}/areas.csv", index = False)

# Generate plots to check nuclei segmentation    
for sample in samples: 
  
  # Define sample subfolder and create folder for plots
  result_folder = f"{output_folder}{sample}"
  if os.path.exists(f"{result_folder}/plots") == False:
    os.mkdir(f"{result_folder}/plots")
    
  # Load filtered anndata object for coordinate info
  adata = sc.read_h5ad(f"{result_folder}/filtered_adata.h5ad")
    
  # Get minimum and maximum x/y coordinates and add padding in each direction
  min_x = int(np.min(adata.obs['pxl_col_in_fullres'])) - 10
  max_x = int(np.max(adata.obs['pxl_col_in_fullres'])) + 10
  min_y = int(np.min(adata.obs['pxl_row_in_fullres'])) - 10
  max_y = int(np.max(adata.obs['pxl_row_in_fullres'])) + 10
    
  # Load DAPI image used for segmentation and normalize 
  img = imread(f'{img_folder}{sample}_DAPI_rollingball.tif')
  img = normalize(img)
    
  # Load nuclei labels 
  labels = np.load(f"{result_folder}/labels.npy")
    
  # Subset image and labels for window containing gene expression
  img = img[min_y:max_y, min_x:max_x]
  labels = labels[min_y:max_y, min_x:max_x]
    
  # Generate lower left coordinates for plot windows 
  y_start = np.random.choice(int(math.floor(np.shape(img)[0]/1000)), size = int(math.floor(np.shape(img)[0]/1000)), replace = False)
  x_start = np.random.choice(int(math.floor(np.shape(img)[1]/1000)), size = int(math.floor(np.shape(img)[1]/1000)), replace = False)
    
  # Initialize plot for 10 windows
  fig, axs = plt.subplots(nrows = 2, ncols = 10, figsize = (100, 20))
  
  # Plot nuclei in random windows  
  for i in range(10):
    min_1 = y_start[i]*1000
    max_1 = min_1 + 1000
    min_2 = x_start[i]*1000
    max_2 = min_2 + 1000
    
    # Shift the minimum x coordinate until at least 50 nuclei are present in the window
    j = i + 1
    while len(np.unique(labels[min_1:max_1,min_2:max_2])) < 50:
      min_2 = x_start[j]*1000
      max_2 = min_2 + 1000
      j += 1
    
    # Plot nuclei    
    axs[0,i].imshow(img[min_1:max_1,min_2:max_2], cmap = "gray")
    axs[0,i].set_title(f"lower left x = {min_2}, y = {min_1}", fontdict={'fontsize': 20})
    axs[0,i].axis("off")
    axs[1,i].imshow(render_label(labels[min_1:max_1,min_2:max_2], img[min_1:max_1,min_2:max_2]))
    axs[1,i].axis("off")
    
    # Go to next window
    i += 1
    
  # Save figure
  fig.savefig(f"{result_folder}/plots/segmentation_windows.pdf", format = "pdf", bbox_inches = "tight")





