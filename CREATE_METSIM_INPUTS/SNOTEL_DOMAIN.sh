############################################################
### Python MetSim State file creator					 ###
### this script opens up snotel files finds the desired  ###
### time period writes netCDF FILES for state inputs     ### 
### for a meteorological simulator						 ###	
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

#######################################
# iterate through months
#######################################
tf=1			#number of time periods to simulate 
for o in range(tf):

	#######################################
	# read in the names of files to format
	# read in the names of .nc output
	#######################################
	with open ("LIST_STATION.csv") as myfile:
		datafile = myfile.read().split('\n')
	
	with open ("LIST_DOMAIN.csv") as myfile:
		outputfile = myfile.read().split('\n')
	
	#######################################
	# define some time variables
	# Used for searching snoTEL data
	# Used 90 days before simulation period
	#######################################
	telstart='2017-01-01'
	telend='2018-09-01'
 
	#######################################
	# IMOPORT .csv files and find time .
	#######################################
	tt=len(datafile)
# 	print(tt)
# 	tt=1
	for z in range(tt):
		data = pd.read_csv(datafile[z], sep=';')

		length=len(data)
		for i in range(length): 
			val = data['date'].values[i]
			coms = (val==telstart)				
			comf = (val==telend)				
			if coms:							#find start of time-series
				start=i
			if comf:							#find end of time-series
				finish=i;

		#######################################
		# extract data from .csv
		#######################################
		length_series=finish-start

		time1 = data['date'].values[start:finish]
		sid = data['site_id'].values[:]
		lon = data['longitude'].values[1]
		lat = data['latitude'].values[1]
		elev = data['elev'].values[1]
		swe = data['snow_water_equivalent'].values[start:finish]
		t_min = data['temperature_min'].values[start:finish]
		t_max = data['temperature_max'].values[start:finish]
		prec = data['precipitation'].values[start:finish]

		#######################################
		# netCDF creation
		#######################################
		outfile=outputfile[z]
# 		outfile+='.'
# 		outfile+=str(telstart)
# 		outfile+='.nc'
		
		ncid = nc4.Dataset(outfile, "w", format="NETCDF4")

		dimid_lon = ncid.createDimension('lon',1)
		dimid_lat = ncid.createDimension('lat',1)
		dimid_m = ncid.createDimension('month',12)

		####################################### Variable: Month
		month_varid = ncid.createVariable('month','i4',('month',))
		length_month=12
		month = [i for i in range(length_month)]

		# Attributes
		# Write data
		month_varid[:] = month

		####################################### Variable: Mask
		mask_varid = ncid.createVariable('mask','i4',('lat','lon',))

		# Attributes
		mask_varid.long_name     = 'domain mask'
		mask_varid.comment       = '0 indicates cell is not active'

		# Write data
		mask_varid[:] = 1

		####################################### Variable: Frac
		frac_varid = ncid.createVariable('frac','d',('lat','lon',))

		# Attributes
		frac_varid.units         = '1'
		frac_varid.long_name     = 'fraction of grid cell that is active'
		frac_varid.FillValue     = 'NaN'

		# Write data
		frac_varid[:] = 1

		####################################### Variable: Elev
		elev_varid = ncid.createVariable('elev','d',('lat','lon',))

		# Attributes
		elev_varid.units         = 'm'
		elev_varid.long_name     = 'gridcell_elevation'
		elev_varid.FillValue    = 'NaN'

		# Write data
		elev_varid[:] = elev
	
		####################################### Variable: Area
		area_varid = ncid.createVariable('area','d',('lat','lon',))

		# Attributes
		area_varid.units         = 'm2'
		area_varid.long_name     = 'area of grid cell'
		area_varid.standardname  = 'area'
		area_varid.FillValue    = 'NaN'

		# Write data
		area_varid[:] = 1

		####################################### Variables: Latitude & Longitude
		lon_varid = ncid.createVariable('lon','d',('lon',))
		lat_varid = ncid.createVariable('lat','d',('lat',))
	
		# Attributes
		lat_varid.names          = '_FillValue'
		lon_varid.names          = '_FillValue'
		lat_varid.value          = 'NaN'
		lon_varid.value          = 'NaN'
	
		# Write data
		lat_varid[:] = lat
		lon_varid[:] = lon

		####################################### Header 
		ncid.License     = 'The file was created by C.Hart, https://github.com/ChristianHart2019'
		ncid.history     = 'Created ' + time.ctime(time.time())

		#####################
		ncid.close()


