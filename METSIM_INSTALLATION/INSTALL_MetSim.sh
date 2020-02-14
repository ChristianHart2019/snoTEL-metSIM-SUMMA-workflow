############################################################
### Python MetSim INSALLATION on Plato					 ###
############################################################

module load python/3.7.4
module load scipy-stack/2019a
module load netcdf
module load hdf5

# DOWNLOAD: pip from link
curl https://bootstrap.pypa.io/get-pip.py -o get-pip.y

pip install scipy --user
pip install netcdf4 --user
pip install numba --user
pip install pandas --user
pip install xarray --user

# DOWNLOAD: MetSim from Github

cd MetSim
python3 setup.py install --user
