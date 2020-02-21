############################################################
### This script creates an initial conditions  files for ###
### SUMMA that is specific to snoTEL site locations	 	 ###	
############################################################

# load all the necessary modules and packages
import pandas as pd
import numpy as np
import re as re
import datetime
import glob
import scipy
import os
import netCDF4 as nc4
import time
import csv
import os
from datetime import date, timedelta
from scipy.io import netcdf
from pathlib import Path

# Define the path with Path() because that avoids issues with how different OS' handle '/' and '\'
list_path = Path('/Users/cjh458/Desktop/LISTS')								# A directory containing list control variables
data_path = Path('/Users/cjh458/Desktop/')									# A path that has an existing initial conditions file
output_path = Path('/Users/cjh458/Desktop/CREATE_SUMMA_INITIALCONDITIONS')	# A path where the output file will be written

# read in a .csv containing snoTEL information and parse into variables
os.chdir(list_path)															# Change directory into the folder containing list data

datafile='LIST_SNOTEL_INFO.csv'												# Open a file that contains the snoTEl information to be written into the attributes file
data = pd.read_csv(datafile, sep=';')										# Read this file using pandas read_csv function to separate variables 
latitude = data['lat'].values[:]											# Extract an Latitude variable from this file
longitude = data['lon'].values[:]											# Extract an Longitude variable from this file
elevation = data['elev'].values[:]											# Extract an Elevation variable from this file
snotelid = data['snotel'].values[:]											# Extract an snotel id variable from this file
hruid = data['hru'].values[:]												# Extract an hru id variable from this file

with open ("LIST_METSIM.csv") as myfile:									# Open list containing names of Metsim-snoTEL simulations
	datafile = myfile.read().split('\n')									# Split up original .CSV into rows

data_len=len(datafile)														# Calculate length of forcing data files (number of rows)
os.chdir(data_path)															# Change directory into the folder that will contain the attribute files

#######################################
# netCDF reading and creation
#######################################

outfile='snotel_initialConditions.nc'										# Create name for an output file
infile='coldState.nc'														# Create name for a data input file

os.chdir(data_path)															# Change directory into the folder that will contain the attribute files
icid = nc4.Dataset(infile, "r", format="NETCDF4")							# Open a file for reading original attributes data
os.chdir(output_path)														# Change directory into the folder that will contain the attribute files
ncid = nc4.Dataset(outfile, "w", format="NETCDF4")							# Open a new file for writing altered attributes data

# Declare dimensions for the new output file				
dimid_hru = ncid.createDimension('hru',data_len)							# Declare hydrological response unit dimension of 631
dimid_gru = ncid.createDimension('gru',data_len)							# Declare geographical response unit dimension of 631
dimid_midSoil = ncid.createDimension('midSoil',8)							# Declare Time-varying variables and parameters at the mid-point of each soil layer
dimid_midToto = ncid.createDimension('midToto',8)							# Declare Time-varying variables and parameters at the mid-point of each layer in the combined soil and snow profile
dimid_midToto = ncid.createDimension('ifcToto',9)							# Declare Time-varying variables and parameters at the interfaces between all layers in the combined soil and snow profile (including top and bottom)
dimid_scalarv = ncid.createDimension('scalarv',1)							# Declare Scalar variables and parameters (degenerate dimension)
dimid_spectral = ncid.createDimension('spectral',1)							# Declare Variables and parameters that vary for different spectral regimes

#######################################
# Declare an hruId value = declare
#######################################
							
hruId_varid = ncid.createVariable('hruId','i4',('hru')) 					# Create hruId variable in the new datafile

# Declare attributes for the new variable 												
hruId_varid.long_name      = 'Ids for hydrological response units'
hruId_varid.units          = '-'
hruId_varid.longname      = 'hru id'

# Assign variable  value
for z in range(data_len):													# loop through the snotel sites	
	hruId_varid[z]= z+1 													# Assign variable a value from the snoTEL information .csv

#######################################
# Declare dt_init 
#######################################
							
dt_init_varid = ncid.createVariable('dt_init','d',('scalarv', 'hru')) 		# Create dt_init variable in the new datafile

raw_data = icid.variables['dt_init']										# Open variable from the "coldState" initial conditions data file for reading 
data = np.copy(raw_data)

# Assign variable  value
for z in range(data_len):													# loop through the snotel sites	
	dt_init_varid[:,z]= data[:,hruid[z]]	 												# Assign variable a value from the coldState information .csv

#######################################
# Declare number of soil layers 
#######################################
							
nSoil_varid = ncid.createVariable('nSoil','i4',('scalarv', 'hru')) 			# Create nSoil variable in the new datafile

raw_data = icid.variables['nSoil']											# Open variable from the "coldState" initial conditions data file for reading 
data = np.copy(raw_data)	

# Assign variable  value
for z in range(data_len):													# loop through the snotel sites	
	nSoil_varid[:,z]= data[:,hruid[z]]														# Assign variable a value from the coldState information .csv

#######################################
# Declare the depth of each layer
#######################################
							
mLayerDepth_varid = ncid.createVariable('mLayerDepth','d',('midToto', 'hru')) # Create mLayerDepth variable in the new datafile

raw_data = icid.variables['mLayerDepth']									# Open variable from the "coldState" initial conditions data file for reading 
data = np.copy(raw_data)													# Make a numpy copy of this variable

# Assign variable  value
for z in range(data_len):													# loop through the snotel sites	
	mLayerDepth_varid[:,z]=data[:,hruid[z]]									# Assign variable a value from the coldState information .csv


######################################
# Declare Height of the layer interface;
#######################################
							
iLayerHeight_varid = ncid.createVariable('iLayerHeight','d',('ifcToto', 'hru')) # Create mLayerDepth variable in the new datafile

raw_data = icid.variables['iLayerHeight']									# Open variable from the "coldState" initial conditions data file for reading 
data = np.copy(raw_data)													# Make a numpy copy of this variable

# Assign variable  value
for z in range(data_len):													# loop through the snotel sites	
	iLayerHeight_varid[:,z]=data[:,hruid[z]]								# Assign variable a value from the coldState information .csv

#######################################
# Declare number of snow layers 
#######################################
							
nSnow_varid = ncid.createVariable('nSnow','i4',('scalarv', 'hru')) 			# Create nSoil variable in the new datafile

raw_data = icid.variables['nSnow']											# Open variable from the "coldState" initial conditions data file for reading 
data = np.copy(raw_data)	

# Assign variable  value
for z in range(data_len):													# loop through the snotel sites	
	nSnow_varid[:,z]=data[:,hruid[z]]  										# Assign variable a value from the coldState information .csv

######################################
# Declare Snow water equivalent;
#######################################
							
scalarSWE_varid = ncid.createVariable('scalarSWE','d',('scalarv', 'hru')) # Create Snow water equivalent variable in the new datafile

raw_data = icid.variables['scalarSWE']										# Open variable from the "coldState" initial conditions data file for reading 
data = np.copy(raw_data)													# Make a numpy copy of this variable

# Assign variable  value
for z in range(data_len):													# loop through the snotel sites	
	scalarSWE_varid[:,z]=data[:,hruid[z]]									# Assign variable a value from the coldState information .csv
	
######################################
# Declare Snow Depth;
#######################################
							
scalarSnowDepth_varid = ncid.createVariable('scalarSnowDepth','d',('scalarv', 'hru')) # Create Snow depth in the new datafile

raw_data = icid.variables['scalarSnowDepth']								# Open variable from the "coldState" initial conditions data file for reading 
data = np.copy(raw_data)													# Make a numpy copy of this variable

# Assign variable  value
for z in range(data_len):													# loop through the snotel sites	
	scalarSnowDepth_varid[:,z]=data[:,hruid[z]]								# Assign variable a value from the coldState information .csv
	
######################################
# Declare Ponded water caused by melt of the "snow without a layer"
#######################################
							
scalarSfcMeltPond_varid = ncid.createVariable('scalarSfcMeltPond','d',('scalarv', 'hru')) # Create Snow Ponded water caused by melt of the "snow without a layer" variable in the new datafile

raw_data = icid.variables['scalarSfcMeltPond']								# Open variable from the "coldState" initial conditions data file for reading 
data = np.copy(raw_data)													# Make a numpy copy of this variable

# Assign variable  value
for z in range(data_len):													# loop through the snotel sites	
	scalarSfcMeltPond_varid[:,z]=data[:,hruid[z]]							# Assign variable a value from the coldState information .csv
	
######################################
# Declare Snow albedo for the entire spectral band
#######################################
							
scalarSnowAlbedo_varid = ncid.createVariable('scalarSnowAlbedo','d',('scalarv', 'hru'))  # Create Snow albedo for the entire spectral band variable in the new datafile

raw_data = icid.variables['scalarSnowAlbedo']								# Open variable from the "coldState" initial conditions data file for reading 
data = np.copy(raw_data)													# Make a numpy copy of this variable

# Assign variable  value
for z in range(data_len):													# loop through the snotel sites	
	scalarSnowAlbedo_varid[:,z]=data[:,hruid[z]]							# Assign variable a value from the coldState information .csv

######################################
# Declare Temperature of the canopy air space;
#######################################
							
scalarCanairTemp_varid = ncid.createVariable('scalarCanairTemp','d',('scalarv', 'hru'))  # Create Temperature of the canopy air space variable in the new datafile

raw_data = icid.variables['scalarCanairTemp']								# Open variable from the "coldState" initial conditions data file for reading 
data = np.copy(raw_data)													# Make a numpy copy of this variable

# Assign variable  value
for z in range(data_len):													# loop through the snotel sites	
	scalarCanairTemp_varid[:,z]=data[:,hruid[z]]							# Assign variable a value from the coldState information .csv
		
######################################
# Temperature of the vegetation canopy
#######################################
							
scalarCanopyTemp_varid = ncid.createVariable('scalarCanopyTemp','d',('scalarv', 'hru')) # Create Temperature of the vegetation canopy variable in the new datafile

raw_data = icid.variables['scalarCanopyTemp']								# Open variable from the "coldState" initial conditions data file for reading 
data = np.copy(raw_data)													# Make a numpy copy of this variable

# Assign variable  value
for z in range(data_len):													# loop through the snotel sites	
	scalarCanopyTemp_varid[:,z]=data[:,hruid[z]]							# Assign variable a value from the coldState information .csv
	
######################################
# Temperature of each layer
#######################################
							
mLayerTemp_varid = ncid.createVariable('mLayerTemp','d',('midToto', 'hru')) # Create Temperature of each layer variable in the new datafile

raw_data = icid.variables['mLayerTemp']										# Open variable from the "coldState" initial conditions data file for reading 
data = np.copy(raw_data)													# Make a numpy copy of this variable

# Assign variable  value
for z in range(data_len):													# loop through the snotel sites	
	mLayerTemp_varid[:,z]=data[:,hruid[z]]									# Assign variable a value from the coldState information .csv

######################################
# Declare Mass of liquid water on the vegetation canopy;
#######################################
							
scalarCanopyLiq_varid = ncid.createVariable('scalarCanopyLiq','d',('scalarv', 'hru'))   # Create Mass of liquid water on the vegetation canopy variable in the new datafile

raw_data = icid.variables['scalarCanopyLiq']								# Open variable from the "coldState" initial conditions data file for reading 
data = np.copy(raw_data)													# Make a numpy copy of this variable

# Assign variable  value
for z in range(data_len):													# loop through the snotel sites	
	scalarCanopyLiq_varid[:,z]=data[:,hruid[z]]								# Assign variable a value from the coldState information .csv
	
######################################
# Matric head of water in the soil;
#######################################
							
mLayerMatricHead_varid = ncid.createVariable('mLayerMatricHead','d',('midSoil', 'hru'))   # Create Matric head of water in the soil variable in the new datafile

raw_data = icid.variables['mLayerMatricHead']								# Open variable from the "coldState" initial conditions data file for reading 
data = np.copy(raw_data)													# Make a numpy copy of this variable

# Assign variable  value
for z in range(data_len):													# loop through the snotel sites	
	mLayerMatricHead_varid[:,z]=data[:,hruid[z]]							# Assign variable a value from the coldState information .csv
	
######################################
# Volumetric fraction of liquid water in each layer
#######################################
							
mLayerVolFracLiq_varid = ncid.createVariable('mLayerVolFracLiq','d',('midToto', 'hru')) # Create Volumetric fraction of liquid water in each layer variable in the new datafile

raw_data = icid.variables['mLayerVolFracLiq']								# Open variable from the "coldState" initial conditions data file for reading 
data = np.copy(raw_data)													# Make a numpy copy of this variable

# Assign variable  value
for z in range(data_len):													# loop through the snotel sites	
	mLayerVolFracLiq_varid[:,z]=data[:,hruid[z]]							# Assign variable a value from the coldState information .csv

######################################
# Declare Mass of ice on the vegetation canopy
#######################################
							
scalarCanopyIce_varid = ncid.createVariable('scalarCanopyIce','d',('scalarv', 'hru'))   # Create Mass of ice on the vegetation canopy variable in the new datafile

raw_data = icid.variables['scalarCanopyIce']								# Open variable from the "coldState" initial conditions data file for reading 
data = np.copy(raw_data)													# Make a numpy copy of this variable

# Assign variable  value
for z in range(data_len):													# loop through the snotel sites	
	scalarCanopyIce_varid[:,z]=data[:,hruid[z]]								# Assign variable a value from the coldState information .csv
	
######################################
# Volumetric fraction of ice in each layer
#######################################
							
mLayerVolFracIce_varid = ncid.createVariable('mLayerVolFracIce','d',('midToto', 'hru')) # Create Volumetric fraction of ice in each layer variable in the new datafile

raw_data = icid.variables['mLayerVolFracIce']								# Open variable from the "coldState" initial conditions data file for reading 
data = np.copy(raw_data)													# Make a numpy copy of this variable

# Assign variable  value
for z in range(data_len):													# loop through the snotel sites	
	mLayerVolFracIce_varid[:,z]=data[:,hruid[z]]							# Assign variable a value from the coldState information .csv

######################################
# Relative aquifer storage -- above bottom of the soil profile
#######################################
							
scalarAquiferStorage_varid = ncid.createVariable('scalarAquiferStorage','d',('scalarv', 'hru'))   # Create Relative aquifer storage -- above bottom of the soil profile variable in the new datafile

raw_data = icid.variables['scalarAquiferStorage']							# Open variable from the "coldState" initial conditions data file for reading 
data = np.copy(raw_data)													# Make a numpy copy of this variable

# Assign variable  value
for z in range(data_len):													# loop through the snotel sites	
	scalarAquiferStorage_varid[:,z]=data[:,hruid[z]]						# Assign variable a value from the coldState information .csv
			
ncid.close()																# close New Datafile
icid.close()																# close "coldState" initial conditions file	




