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
# Summary: Identify hippocampal nuclei
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
output_dir = "/path/to/hippocampus/output/"
roi_dir = "/path/to/rois/"

# Select sample (B9, A9, A34995E) and create subfolder
sample = "B9"
if os.path.exists(f"{output_dir}{sample}") == False:
  os.mkdir(f"{output_dir}{sample}")
  
# Load ROI 
roi = read_roi_file(f"{roi_dir}hippocampus-{sample}.roi")
roi = roi[f"hippocampus-{sample}"]

# Extract x and y coordinates from the ROI 
x_coords = roi['x']
y_coords = roi['y']
coordinates = list(zip(x_coords, y_coords))

# Create a Polygon
polygon = Polygon(coordinates)

# Load GDF from nuclei segmentation 
gdf_nuclei = gpd.read_file(f"{input_dir}{sample}/gdf/gdf.shp")

# Identify nuclei within the ROI
in_hippocampus = gdf_nuclei[gdf_nuclei.geometry.within(polygon)]

# Save hippocampal nuclei GDF
if os.path.exists(f"{output_dir}{sample}/gdf_hippocampus") == False:
  os.mkdir(f"{output_dir}{sample}/gdf_hippocampus")
in_hippocampus.to_file(f"{output_dir}{sample}/gdf_hippocampus/gdf.shp")

# Save filtered nuclei to CSV
in_hippocampus[['id']].to_csv(f"{output_dir}{sample}/nuclei_in_hippocampus.csv", index=False)









