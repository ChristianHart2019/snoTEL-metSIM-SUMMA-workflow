#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=0-01:00

module load python/3.7.4
module load scipy-stack/2019a
module load netcdf
module laod hdf5

python3 windspd_NLDAS3.sh
    
