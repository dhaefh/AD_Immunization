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
# Summary: Calculate distance to nearest parenchymal amyloid for each nucleus
#
#-------------------------------------------------------------------------------
# Initialization

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
from stardist import random_label_cmap
from stardist import render_label
from shapely.geometry import Polygon, Point
from scipy import sparse
from matplotlib.colors import ListedColormap
import os
import math
import pickle
from shapely.ops import nearest_points
import matplotlib.patches as patches
from read_roi import read_roi_zip
from read_roi import read_roi_file

# Define input and output folders 
input_dir = "/path/to/preprocessing/output/folder/"
output_dir = "/path/to/parenchymal/amyloid/distance/output/"
roi_dir = "/path/to/rois/"

# Select sample (B9, A9, A34995E) and create subfolder
sample = "B9"
if os.path.exists(f"{output_dir}{sample}") == False:
  os.mkdir(f"{output_dir}{sample}")
  
# Load ROI zip file
rois = read_roi_zip(f"{roi_dir}parenchymal-amyloid-{sample}.zip")

# Convert odict_values object to list in order to iterate through ROIs
roi_list = list(rois.values())

# Create a list of Polygon geometries
roi_geometries = []
for i in range(len(roi_list)): 
  
  # Extract current ROI 
  roi = roi_list[i]
    
  # Extract coordinates from ROI
  x_coords = roi['x']
  y_coords = roi['y']
  coordinates = list(zip(x_coords, y_coords))
    
  # Create a Polygon
  polygon = Polygon(coordinates)
    
  # Fix invalid polygons
  if not polygon.is_valid:
    print("Fixing self-intersection...")
    polygon = polygon.buffer(0)
  if not polygon.is_valid:
    print("Polygon still not valid")

  roi_geometries.append(polygon)

# Create a GeoDataFrame
gdf_roi = gpd.GeoDataFrame(geometry=roi_geometries)
gdf_roi['id'] = [f"roi_{i+1}" for i in range(len(gdf_roi))]

# Save ROI GDF
if os.path.exists(f"{output_dir}{sample}/gdf_roi") == False:
  os.mkdir(f"{output_dir}{sample}/gdf_roi")
gdf_roi.to_file(f"{output_dir}{sample}/gdf_roi/gdf.shp")

# Calculate union of multiple geometries (MultiPolygon object) 
roi_union = gdf_roi['geometry'].unary_union

# Load GDF from nuclei segmentation 
gdf_nuclei = gpd.read_file(f"{input_dir}{sample}/gdf/gdf.shp")

# Define function to calculate shortest distance to ROI for each nucleus
def nearest_distance(gdf):
  # Identify the point in the ROI closest to the nucleus (returns the nucleus itself if it is within the bounds of an ROI)
  nearest_geom = nearest_points(gdf['geometry'], roi_union)[1] 
    
  # Calculate distance from nucleus to closest point in ROI (zero if nucleus is within ROI)
  return gdf['geometry'].distance(nearest_geom)

# Calculate distance to ROI for each nucleus
gdf_nuclei['roi_dist'] = gdf_nuclei.apply(nearest_distance, axis=1)

# Save nuclei GDF with distance to parenchymal amyloid
if os.path.exists(f"{output_dir}{sample}/gdf_roi_distance") == False:
  os.mkdir(f"{output_dir}{sample}/gdf_roi_distance")
gdf_nuclei.to_file(f"{output_dir}{sample}/gdf_roi_distance/gdf.shp")

# Save distances to csv
dist = pd.DataFrame(gdf_nuclei[['id','roi_dist']])
dist.to_csv(f"{output_dir}{sample}/distance.csv", index = False)





