############################################################
### Python concatenates windspeed data from nldas		 ###
### and attempts to create time series files for snoTEl  ###
### that will be used as an input for SUMMA				 ###	
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

# os.chdir('/Users/cjh458/Desktop/LISTS')	# Change directory into the folder containing lists of names
os.chdir('/datastore/GLOBALWATER/CommonData/NLDAS_SUMMA/SNOTEL_NLDAS_WDSP')	# Change directory into the folder containing lists of names
	
#######################################
# iterate through time periods
#######################################
with open ("LIST_HRU.csv") as myfile:								# Declare SNOTEL file names and corresponding HRU # to find
		hrufile = myfile.read().split('\n')							# Split up original .CSV into rows

with open ("LIST_HRU_SNOTEL.csv") as myfile:						# Declare SNOTEL file names and corresponding HRU # to find
		snotelfile = myfile.read().split('\n')						# Split up original .CSV into rows

snotellength=5#len(hrufile)											# Calculate length of these data files
for o in range(snotellength):										# Iterate through files

	print(o)	
	
	time_data=[]													# create time data variable
	wind_data=[]													# create wind data variable

	#######################################
	# read in the names of metsim output files
	# read in the names of SUMMA .nc output files
	#######################################
# 	os.chdir('/Users/cjh458/Desktop/LISTS')							# Change directory into the folder containing lists of names
	os.chdir('/datastore/GLOBALWATER/CommonData/NLDAS_SUMMA/SNOTEL_NLDAS_WDSP')	 
	with open ("LIST_NLDAS.csv") as myfile:							# Open list containing names of Metsim-snoTEL simulations
		nldasfile = myfile.read().split('\n')						# Split up original .CSV into rows
	os.chdir('/datastore/GLOBALWATER/CommonData/NLDAS_SUMMA/')		# Change directory into the folder containing NLDAS data
# 	os.chdir('/Users/cjh458/Desktop/NLDAS')							# Change directory into the folder containing NLDAS data		
	
	#######################################
	# open NLDAS .nc files 
	#######################################
	nldaslength=len(nldasfile)	# 	tt=1							# Calculate number of NLDAS data files (number of rows)
	for z in range(nldaslength):									# Iterate through files
		
		#######################################
		# netCDF opening and creation
		#######################################		
# 		mcid = nc4.Dataset(nldasfile[z], "r", format="NETCDF4")		
		filename=nldasfile[z]										# Assign a file name variable to a FILE name from the nldas name list file
# 		filename='nldasForcing_201701.nc'
		mcid = nc4.Dataset(filename, "r", format="NETCDF4")			# Open the NLDAS data file for reading	
			
		#######################################
		# Format windspeed @ time stamp 
		#######################################		
		raw_data = mcid.variables['time']							# Assign a raw_data variable to the time variable from the NLDAS data
		time_data = np.append(time_data,raw_data)					# Append the raw_data to the entire time time-series
		raw_data = mcid.variables['windspd']						# Assign a raw_data variable to the windspd variable from the NLDAS data
		wind_data = np.append(wind_data,raw_data[:,hrufile[o]])		# Append the raw_data to the entire windspd time-series
		
# 		print(len(hrufile))
		mcid.close()												# close NLDAS Datafile

# 	os.chdir('/Users/cjh458/Desktop/NLDAS')							# Change directory into the folder containing lists of names
	os.chdir('/datastore/GLOBALWATER/CommonData/NLDAS_SUMMA')	 	# Change directory into the folder containing NLDAS data

	#######################################
	# netCDF creating a new file
	#######################################	
	outfile=snotelfile[o]										    # Create name for an output file from and the existing snoTEL file
	outfile+='.wind.nc'												# Concatenate a wind variable identifier to the original file name
	ncid = nc4.Dataset(outfile, "w", format="NETCDF4")				# Create a file for writing NLDAS data using this new name
	# Declare dimensions for the new output file				
	dimid_T = ncid.createDimension('time',14592)
	dimid_hru = ncid.createDimension('hru',1)
	# Declare attributes for the new variable 				
	wind_varid = ncid.createVariable('windspd','f',('time','hru'))	# Create variable used to assign the Windspd time-series to
	wind_varid.long_name      = 'wind speed at the measurement height' # Create long_name attribute for the windspd variable
	wind_varid.units          = 'm s-1'								# Create units attribute for the windspd variable
	wind_varid.FillValue      = '-999'								# Create a FillValue attribute for the windspd variable
	# Assign variable  value
	wind_varid[:]=wind_data											# Assign NLDAS data to file
	# Declare attributes for the new variable 				
	time_varid = ncid.createVariable('time','d',('time','hru'))		# Create variable used to assign the time time-series to
	time_varid.long_name      = 'Observation time'					# Create long_name attribute for the time variable
	time_varid.units          = 'days since 1990-01-01 00:00:00'	# Create units attribute for the time variable
	time_varid.FillValue      = '-999'								# Create a FillValue attribute for the time variable
	# Assign variable a dummy value
	time_varid[:]=time_data											# Assign NLDAS data to file

	ncid.close()													# Close WINDSPD Datafile	
