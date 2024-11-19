#!/bin/bash
#SBATCH --account p31535
#SBATCH --partition normal
#SBATCH --job-name ptau-rollingball
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 32G
#SBATCH --time 6:00:00
#SBATCH --output /projects/p31535/alex/fiji-image-analysis/an1792-rebuttal/%x.oe%j.log
#SBATCH --verbose

module load fiji

ImageJ-linux64 --headless -macro /projects/b1169/alex/scripts/AN1792-rebuttal/hippocampal/background-subtraction-ptau.ijm