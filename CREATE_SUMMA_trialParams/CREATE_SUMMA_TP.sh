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
output_path = Path('/Users/cjh458/Desktop/CREATE_SUMMA_trialParams')		# A path where the output file will be written

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

outfile='snotel_trialParams.nc'										# Create name for an output file
infile='trialParams.nc'														# Create name for a data input file

os.chdir(data_path)															# Change directory into the folder that will contain the attribute files
icid = nc4.Dataset(infile, "r", format="NETCDF4")							# Open a file for reading original attributes data
os.chdir(output_path)														# Change directory into the folder that will contain the attribute files
ncid = nc4.Dataset(outfile, "w", format="NETCDF4")							# Open a new file for writing altered attributes data

# Declare dimensions for the new output file				
dimid_hru = ncid.createDimension('hru',data_len)							# Declare hydrological response unit dimension of 631
dimid_gru = ncid.createDimension('gru',data_len)							# Declare geographical response unit dimension of 631

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
# Declare an gruId value = declare
#######################################
							
gruId_varid = ncid.createVariable('gruId','i4',('gru')) 					# Create gruId variable in the new datafile

# Declare attributes for the new variable 												
gruId_varid.long_name      = 'Ids for geographical response units'
gruId_varid.units          = '-'
gruId_varid.longname      = 'gru id'

# Assign variable  value
for z in range(data_len):													# loop through the snotel sites	
	gruId_varid[z]= z+1 


#######################################
# Declare an critSoilTranspire value = declare
#######################################
							
critSoilTranspire_varid = ncid.createVariable('critSoilTranspire','d',('hru')) 					# Create critSoilTranspire variable in the new datafile

# Assign variable  value
for z in range(data_len):													# loop through the snotel sites	
	critSoilTranspire_varid[z]= 0.175


#######################################
# Declare an theta_res value = declare
#######################################
							
theta_res_varid = ncid.createVariable('theta_res','d',('hru')) 					# Create critSoilTranspire variable in the new datafile

# Assign variable  value
for z in range(data_len):													# loop through the snotel sites	
	theta_res_varid[z]= 0.139


#######################################
# Declare an theta_sat value = declare
#######################################
							
theta_sat_varid = ncid.createVariable('theta_sat','d',('hru')) 					# Create critSoilTranspire variable in the new datafile

# Assign variable  value
for z in range(data_len):													# loop through the snotel sites	
	theta_sat_varid[z]= 0.55


ncid.close()																# close New Datafile
icid.close()																# close "coldState" initial conditions file	




