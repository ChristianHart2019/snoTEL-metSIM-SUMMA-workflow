############################################################
### This script creates attribute files for SUMMA that   ###
### are specific to snoTEL site locations				 ###	
############################################################

# load all the necessary modules and packages
import pandas as pd
import numpy as np
import re as re
import datetime
from datetime import date, timedelta
import glob
import scipy
import os
from scipy.io import netcdf
import netCDF4 as nc4
import time
import csv
import os


#######################################
# read in a .csv containing snoTEL
# information and parse into variables
#######################################
os.chdir('/Users/cjh458/Desktop/LISTS')							# Change directory into the folder containing list data

datafile='LIST_SNOTEL_INFO.csv'									# Open a file that contains the snoTEl information to be written into the attributes file
data = pd.read_csv(datafile, sep=';')							# Read this file using pandas read_csv function to separate variables 
latitude = data['lat'].values[:]								# Extract an Latitude variable from this file
longitude = data['lon'].values[:]								# Extract an Longitude variable from this file
elevation = data['elev'].values[:]								# Extract an Elevation variable from this file
snotelid = data['snotel'].values[:]								# Extract an snotel id variable from this file
hruid = data['hru'].values[:]									# Extract an hru id variable from this file

#######################################
# read in the names of metsim output files
#######################################
with open ("LIST_METSIM.csv") as myfile:					# Open list containing names of Metsim-snoTEL simulations
	datafile = myfile.read().split('\n')					# Split up original .CSV into rows

tt=len(datafile)													# Calculate length of forcing data files (number of rows)
																# Iterate through the forcing file names
# 	os.mkdir('/Users/cjh458/Desktop/CREATE_snoTEL_ATTRIBUTES')	# Create directory that will contain the attribute files
os.chdir('/Users/cjh458/Desktop/CREATE_SUMMA_ATTRIBUTES')		# Change directory into the folder that will contain the attribute files

#######################################
# netCDF reading and creation
#######################################

# snotelid[z]='{0:g}'.format(snotelid[z])							# remove trailing zeros from the imported floating point value
# outfile='snotel_' + str(snotelid[z]) + '_20170101-20180831.nc.attributes.nc' # Create name for an output file 
outfile='snotel_attributes.nc'									# Create name for an output file
attribfile='attributes.nc'										# Create file identifier for the original attributes file

acid = nc4.Dataset(attribfile, "r", format="NETCDF4")			# Open a file for reading original attributes data
ncid = nc4.Dataset(outfile, "w", format="NETCDF4")				# Open a new file for writing altered attributes data
# Declare dimensions for the new output file				
dimid_hru = ncid.createDimension('hru',tt)						# Declare hydrological response unit dimension of 1
dimid_gru = ncid.createDimension('gru',tt)						# Declare geographical response unit dimension of 1


#######################################
# Declare an HRUarea value of 1
#######################################
							
HRUarea_varid = ncid.createVariable('HRUarea','d',('hru')) 		# Create time variable in the new datafile

# Declare attributes for the new variable 				
HRUarea_varid.long_name      = 'Area of each HRU'
HRUarea_varid.units          = 'm^2'
HRUarea_varid.FillValue      = '-999'

# Assign variable  value
for z in range(tt):												# loop through the snotel sites
	HRUarea_varid[z]=1											# Assign variable a constant value	


#######################################
# Declare an contourLength a value of 1
#######################################
							
contourLength_varid = ncid.createVariable('contourLength','d',('hru')) 	# Create time variable in the new datafile

# Declare attributes for the new variable 				
contourLength_varid.long_name      = 'ContourLength of HRU'
contourLength_varid.units          = 'm'

# Assign variable  value
for z in range(tt):												# loop through the snotel sites
	contourLength_varid[z]=1									# Assign variable a constant value	


#######################################
# Declare an downHRUindex value of 0
#######################################
							
downHRUindex_varid = ncid.createVariable('downHRUindex','i4',('hru')) 	# Create time variable in the new datafile

# Declare attributes for the new variable 				
downHRUindex_varid.long_name      = 'Id of downslope HRU (0 = basin outlet)'
downHRUindex_varid.units          = '-'

# Assign variable  value
for z in range(tt):												# loop through the snotel sites
	downHRUindex_varid[z]=0											# Assign variable a constant value


#######################################
# Declare an elevation value = file
#######################################
							
elevation_varid = ncid.createVariable('elevation','f',('hru')) 	# Create time variable in the new datafile

# Declare attributes for the new variable 				
elevation_varid.long_name      = 'NLDAS elevation'
elevation_varid.units          = 'meter'
elevation_varid.FillValue      = '-999.f'

# Assign variable  value
for z in range(tt):												# loop through the snotel sites	
	elevation_varid[z]=elevation[z]									# Assign variable a value from the snoTEL information .csv


#######################################
# Declare an hruId value = file
#######################################
							
hruId_varid = ncid.createVariable('hruId','i4',('hru')) 		# Create time variable in the new datafile

# Declare attributes for the new variable 				
raw_data = acid.variables['hru2gruId']							# Open existing hru2gruId variable 		
h2g_data = np.copy(raw_data)									# Copy this variable to a numpy array
hruId_varid.long_name      = 'hru id'
hruId_varid.units          = '-'

# Assign variable  value
for z in range(tt):												# loop through the snotel sites	
	hruId_varid[z]= z+1 #103473 #h2g_data[hruid[z]]#hruid[z]			# Assign variable a value from the snoTEL information .csv


#######################################
# Declare an latitude value of 1
#######################################
							
latitude_varid = ncid.createVariable('latitude','d',('hru')) 	# Create time variable in the new datafile

# Declare attributes for the new variable 				
latitude_varid.long_name      = 'Latitude of HRU\'s centriod point'
latitude_varid.units          = 'decimal degree north'

# Assign variable  value
for z in range(tt):												# loop through the snotel sites	
	latitude_varid[z]=latitude[z]									# Assign variable a value from the snoTEL information .csv


#######################################
# Declare an longitude value of 1
#######################################
							
longitude_varid = ncid.createVariable('longitude','d',('hru')) 	# Create time variable in the new datafile

# Declare attributes for the new variable 				
longitude_varid.long_name      = 'longitude of HRU\'s centriod point'
longitude_varid.units          = 'decimal degree east'

# Assign variable  value
for z in range(tt):												# loop through the snotel sites
	longitude_varid[z]=longitude[z]									# Assign variable a value from the snoTEL information .csv			


#######################################
# Declare an mHeight value of 3
#######################################
							
mHeight_varid = ncid.createVariable('mHeight','d',('hru')) 		# Create time variable in the new datafile

# Declare attributes for the new variable 				
mHeight_varid.long_name      = 'Measurement height above bare ground'
mHeight_varid.units          = 'm'

# Assign variable  value
for z in range(tt):												# loop through the snotel sites
	mHeight_varid[z]=3												# Assign variable a constant value


#######################################
# Declare an slopeTypeIndex = existing
#######################################
							
slopeTypeIndex_varid = ncid.createVariable('slopeTypeIndex','i4',('hru')) 	# Create time variable in the new datafile

# Declare attributes for the new variable 				
slopeTypeIndex_varid.long_name      = 'index defining slope'
slopeTypeIndex_varid.units          = '-'
slopeTypeIndex_varid.FillValue      = '-999'

# Assign variable  value
raw_data = acid.variables['slopeTypeIndex']						# Open existing slopeTypeIndex variable 		
sti_data = np.copy(raw_data)									# Copy this variable to a numpy array
for z in range(tt):												# loop through the snotel sites
	slopeTypeIndex_varid[z] = sti_data[hruid[z]]					# assign new variable existing value for the HRU	


#######################################
# Declare an tan_slope value of 0.1
#######################################
							
tan_slope_varid = ncid.createVariable('tan_slope','d',('hru')) 	# Create time variable in the new datafile

# Declare attributes for the new variable 				
tan_slope_varid.long_name      = 'Average tangent slope of HRU'
tan_slope_varid.units          = 'm m-1'

# Assign variable  value
for z in range(tt):												# loop through the snotel sites
	tan_slope_varid[z]=0.1											# Assign variable a constant value


#######################################
# Declare an vegTypeIndex = existing
#######################################
							
vegTypeIndex_varid = ncid.createVariable('vegTypeIndex','i4',('hru')) 	# Create time variable in the new datafile

# Declare attributes for the new variable 				
vegTypeIndex_varid.long_name      = 'Index defining vegetation type'
vegTypeIndex_varid.units          = '-'
vegTypeIndex_varid.FillValue      = '-999'

# Assign variable  value
raw_data = acid.variables['vegTypeIndex']						# Open existing vegTypeIndex variable 		
vti_data = np.copy(raw_data)									# Copy this variable to a numpy array
for z in range(tt):												# loop through the snotel sites
	vegTypeIndex_varid[z] = vti_data[hruid[z]]						# assign new variable existing value for the HRU


#######################################
# Declare an hru2gruId = existing
#######################################
							
hru2gruId_varid = ncid.createVariable('hru2gruId','i4',('hru')) # Create time variable in the new datafile

# Declare attributes for the new variable 				
hru2gruId_varid.long_name      = 'Id of GRU to which the HRU belongs'
hru2gruId_varid.units          = '-'

# Assign variable  value
for z in range(tt):												# loop through the snotel sites
	hru2gruId_varid[z] = z+1 #103473 #h2g_data[hruid[z]] #hruid[z]				# assign new variable existing value for the HRU


#######################################
# Declare an gruId = existing
#######################################
							
gruId_varid = ncid.createVariable('gruId','i4',('gru')) 		# Create time variable in the new datafile

# Declare attributes for the new variable 				
gruId_varid.long_name      = 'Id of GRU to which the HRU belongs'
gruId_varid.units          = '-'

# Assign variable  value
raw_data = acid.variables['gruId']								# Open existing gruId variable 		
h2g_data = np.copy(raw_data)									# Copy this variable to a numpy array
for z in range(tt):												# loop through the snotel sites
	gruId_varid[z] = z+1 #103473 #h2g_data[hruid[z]]		#hruid[z				# assign new variable existing value for the HRU


#######################################
# Declare an soilTypeIndex
#######################################
							
soilti_varid = ncid.createVariable('soilTypeIndex','i4',('hru')) 		# Create time variable in the new datafile

# Declare attributes for the new variable 				
soilti_varid.long_name      = 'Index defining soil type'
soilti_varid.units          = '-'
soilti_varid.FillValue      = '-999'

# Assign variable  value
raw_data = acid.variables['soilTypeIndex']						# Open existing gruId variable 		
soilti_data = np.copy(raw_data)									# Copy this variable to a numpy array
for z in range(tt):												# loop through the snotel sites
	soilti_varid[z] = soilti_data[hruid[z]]							# assign new variable existing value for the HRU


acid.close()													# close netCDF dataFile
ncid.close()													# close netCDF dataFile
	
	