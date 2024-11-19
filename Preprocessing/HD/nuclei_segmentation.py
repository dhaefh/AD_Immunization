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
# Summary: Nuclei segmentation with StarDist
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

# Define general output folder
output_folder = "/path/to/preprocessing/output/folder/"

# Define path to DAPI images
img_folder = '/path/to/DAPI/images/'

# Define list of samples
samples = ['B9', 'A9', 'A34995E']

# Load the pretrained model
model = StarDist2D.from_pretrained('2D_versatile_fluo')

# Segment nuclei for each sample
for sample in samples:
  
  # Create subfolder
  result_folder = f'{output_folder}{sample}'
  if os.path.exists(result_folder) == False:
    os.mkdir(result_folder)
    
  # Load DAPI image
  img = imread(f'{img_folder}{sample}_DAPI_rollingball.tif')
  
  # Normalize image
  img = normalize(img)

  # Identify nuclei using default probability and NMS thresholds
  labels, polys = model.predict_instances_big(img, axes='YX', block_size=4096, min_overlap=128, normalizer=None, n_tiles=(4,4))
    
  # Save output
  np.save(f"{result_folder}/labels", labels)
  file = open(f"{result_folder}/polys.pkl","wb")
  pickle.dump(polys,file)
  file.close()








