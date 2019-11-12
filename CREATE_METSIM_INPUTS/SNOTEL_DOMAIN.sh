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
# define some time variables
# Used for searching snoTEL data
#######################################
with open ("LIST_START.csv") as myfile:
    startfile = myfile.read().split('\n')	#Read list containing the end of the state data time-series, which is = to start of the simulation

with open ("LIST_END.csv") as myfile:
    endfile = myfile.read().split('\n')	#Read list containing the end of the state data time-series, which is = to start of the simulation

timedata = pd.read_csv('LIST_TIME.csv', sep=',') # read in general information about the start and end of the simulation period

month = timedata['month'].values[:]
startm = timedata['startm'].values[:]
starty = timedata['starty'].values[:]
endy = timedata['endy'].values[:]
endd = timedata['endd'].values[:] 
endm = timedata['endm'].values[:] 

#######################################
# iterate through months
#######################################
tf=13			#number of months to simulate 
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
	telstart=str(startfile[o])
	telend=str(endfile[o])
# 	telstart='2017-08-01'
# 	telend='2017-09-01'
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
		outfile+='.'
		outfile+=str(starty[o])
		outfile+='.'
		outfile+=str(startm[o])
		outfile+='.nc'
		
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
	
		lat1 = [30.0312500000000, 30.0937500000000, 30.1562500000000, 30.2187500000000, 30.2812500000000, 30.3437500000000, 30.4062500000000, 30.4687500000000, 30.5312500000000]	
		lon1 = [-100.031250000000, -99.9687500000000, -99.9062500000000, -99.8437500000000, -99.7812500000000, -99.7187500000000, -99.6562500000000, -99.5937500000000, -99.5312500000000]
	
		length_lon=len(lon1)
		for i in range(length_lon): 
			lat1[i] = lat-(0.001*i)
			lon1[i] = lon-(0.001*i)
	
		lat_varid[:] = lat
		lon_varid[:] = lon


		####################################### Header 
		ncid.License     = 'The file was created by C.Hart, https://github.com/ChristianHart2019'
		ncid.history     = 'Created ' + time.ctime(time.time())

		#####################
		ncid.close()


