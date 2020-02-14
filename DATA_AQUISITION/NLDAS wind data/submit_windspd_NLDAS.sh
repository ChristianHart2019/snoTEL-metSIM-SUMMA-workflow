#!/bin/bash
#SBATCH --nodes=1               
#SBATCH --time=0-01:00

#load require modules
module load python/3.7.4
module load scipy-stack/2019a
module load netcdf
module laod hdf5

# run the win acquistion script
python3 windspd_NLDAS3.sh
    
