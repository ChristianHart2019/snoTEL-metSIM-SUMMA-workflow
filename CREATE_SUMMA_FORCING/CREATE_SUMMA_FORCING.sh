
############################################################
### This script formats output variables from MetSim     ###
### and NLDAS netCDF files for use with SUMMA			 ###	
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

os.chdir('/Users/cjh458/Desktop/LISTS')								# Change directory into the folder containing lists of names
	
#######################################
# iterate through time periods
#######################################
tf=1																# Declare number of time periods to use -- or files to create
for o in range(tf):													# Iterate through these periods
	
	#######################################
	# read in the names of metsim output files
	#######################################
	with open ("LIST_METSIM1.csv") as myfile:						# Open list containing names of Metsim-snoTEL simulations
		datafile = myfile.read().split('\n')						# Split up original .CSV into rows

	#######################################
	# open MetSim and NLDAS .nc files 
	#######################################
	tt=1#len(datafile)	# 	tt=1									# Calculate length of MetSim data (number of rows)
	for z in range(tt):												# Iterate through the MetSim file names
		os.chdir('/Users/cjh458/Desktop/SNOTEL_FORCING_DATA')		# Change directory into the folder with the MetSim output data		
		
		#######################################
		# netCDF opening and creation
		#######################################
		mcid = nc4.Dataset(datafile[z], "r", format="NETCDF4")		# Open the MetSim data file for reading
		outfile=datafile[z]											# Create name for an output file from and the existing input file
		outfile+='.nc'	
		windfile=datafile[z]										# Create name for an output file from and the existing input file
		windfile+='.wind.nc'
		print(windfile)												# Create name for an output file from and the existing input file
		wcid = nc4.Dataset(windfile, "r", format="NETCDF4")			# Open the wind data file for reading
		ncid = nc4.Dataset(outfile, "w", format="NETCDF4")			# Open a new file for writing Metsim and NLDAS data
		# Declare dimensions for the new output file				
		dimid_T = ncid.createDimension('time',len(mcid.variables['time'])) # Declare time dimension that is the same length of the data file
		dimid_hru = ncid.createDimension('hru',1)					# Declare hydrological response unit dimension of 1
		
		#######################################
		# attain and format the windspeed
		#######################################
		raw_data = wcid.variables['windspd']						# Open the time variable 		
		data = np.copy(raw_data)									# Copy this data into numpy format
		wind_varid = ncid.createVariable('windspd','f',('time','hru'))	# Create time variable in the new datafile

		# Declare attributes for the new variable 				
		wind_varid.long_name      = 'wind speed at the measurement height'
		wind_varid.units          = 'm s-1'
		wind_varid.FillValue      = '-999'
		
		# Assign variable  value
		wind_varid[:]=1#data											# Assign NLDAS data to file

		#######################################
		# Format the time stamp 
		#######################################
		raw_data = mcid.variables['time']							# Open the time variable 
		data = np.copy(raw_data)									# Make a copy of this variable
		timedata=((data))											# Convert days since 1990-01-01 00:00:00 to hours since 2000-01-01 00:00:00			
		time_varid = ncid.createVariable('time','d',('time','hru'))	# Create time variable in the new datafile

		# Declare variable attributes
		time_varid.long_name      = 'Observation time'
		time_varid.units          = 'hours since 2000-01-01 00:00:00'
		time_varid.calendar       = 'standard'

		# Write data to the variable
		time_varid[:] = timedata									# Update the original time stamp to the converted one
	
		#######################################
		# Change names and attrib. of radiation 
		#######################################
		raw_data = mcid.variables['longwave']						# Open the longwave variable 
		raw_data1 = mcid.variables['shortwave']						# Open the shortwave variable 
		data = np.copy(raw_data)									# Make a copies of these variables
		data1 = np.copy(raw_data)									# Make a copies of these variables
		
		# Create longwave variable in the new datafile
		LWRadAtm_varid = ncid.createVariable('LWRadAtm','f',('time','hru'))	

		# Declare variable attributes
		LWRadAtm_varid.long_name      = 'downward longwave radiation at the upper boundary'
		LWRadAtm_varid.units          = 'W m-2'
		LWRadAtm_varid.FillValue      = '-999'
		
		# Write data to the variable
		LWRadAtm_varid[:] = data
											
		# Create shortwave variable in the new datafile
		SWRadAtm_varid = ncid.createVariable('SWRadAtm','f',('time','hru'))	
		
		# Declare variable attributes
		SWRadAtm_varid.long_name      = 'downward shortwave radiation at the upper boundary'
		SWRadAtm_varid.units          = 'W m-2'
		SWRadAtm_varid.FillValue      = '-999'
			
		# Write data to the variable						
		SWRadAtm_varid[:] = data1	
		
		#######################################
		# Format the airpres
		#######################################
		raw_data = mcid.variables['air_pressure']					# Open the air_pressure variable 
		data = np.copy(raw_data)									# Make a copy of this variable
		presdata=(data*1000)										# Convert kPa to Pa			
		
		pres_varid = ncid.createVariable('airpres','f',('time','hru'))	# Create airpressure variable in the new datafile

		# Declare variable attributes
		pres_varid.long_name      = 'air pressure at the measurement height'
		pres_varid.units          = 'Pa'
		pres_varid.FillValue      = '-999'

		# Write data to the variable
		pres_varid[:] = presdata									# Update the original air pressure to the converted one
									
		#######################################
		# Format the airtemp
		#######################################
		raw_data = mcid.variables['temp']							# Open the air_pressure variable 
		data = np.copy(raw_data)									# Make a copy of this variable
		tempdata=(data+237.15)										# Convert C to K			
		
		temp_varid = ncid.createVariable('airtemp','f',('time','hru'))	# Create airtemp variable in the new datafile

		# Declare variable attributes
		temp_varid.long_name      = 'air temperature at the measurement height'
		temp_varid.units          = 'K'
		temp_varid.FillValue      = '-999'

		# Write data to the variable
		temp_varid[:] = tempdata									# Update the original air temperature to the converted one
									
		#######################################
		# Create data_step varaiable
		#######################################			
		
		datastep_varid = ncid.createVariable('data_step','d')		# Create data_step variable in the new datafile

		# Declare variable attributes
		datastep_varid.long_name      = 'data step length in seconds'
		datastep_varid.units          = 'seconds'

		# Write data to the variable
		datastep_varid[:] = 3600									# carry over the original data_step to the new one
	
		#######################################
		# Create hruId varaiable
		#######################################			
		
		hruId_varid = ncid.createVariable('hruId','i',('hru'))		# Create hruId variable in the new datafile

		# Declare variable attributes
		hruId_varid.long_name      = 'hru id'
		hruId_varid.units          = '-'

		# Write data to the variable
		hruId_varid[:] = 1											# Update the original hruId to the new one

		#######################################
		# Format Precipitation variable
		#######################################
		raw_data = mcid.variables['prec']							# Open the prec variable 
		data = np.copy(raw_data)									# Make a copy of this variable
		precdata=(data*3600)										# Convert mm/timestep to kg m-2 s-1			
		
		prec_varid = ncid.createVariable('pptrate','d',('time','hru'))	# Create pptrate variable in the new datafile

		# Declare variable attributes
		prec_varid.long_name      = 'Precipitation rate'
		prec_varid.units          = 'kg m-2 s-1'
		prec_varid.FillValue      = '-999'

		# Write data to the variable
		prec_varid[:] = precdata									# Update the original precipitation data and add to the converted dataset 
		
		#######################################
		# Format Spechum variable
		#######################################
		raw_data = mcid.variables['vapor_pressure']					# Open the vapor pressure variable 
		raw_data1 = mcid.variables['air_pressure']					# Open the air pressure variable 
		data = np.copy(raw_data)									# Make a copies of these variables
		data1 = np.copy(raw_data1)									# Make a copies of these variables
											
		spechumdata=(0.622*data)/(data1-data*(1-0.622))				# Convert air and vapor pressure to specific humidity 		

		spechum_varid = ncid.createVariable('spechum','f',('time','hru'))	# Create specific humidity variable in the new datafile

		# Declare variable attributes
		spechum_varid.long_name      = 'specific humidity at the measurement height'
		spechum_varid.units          = 'g g-1'
		spechum_varid.FillValue      = '-999'

		# Write data to the variable
		spechum_varid[:] = spechumdata								# Update the original specific humidity data and add to the converted dataset 
		
		mcid.close()												# close MetSim Datafile
		ncid.close()												# close MetSim Datafile
		
