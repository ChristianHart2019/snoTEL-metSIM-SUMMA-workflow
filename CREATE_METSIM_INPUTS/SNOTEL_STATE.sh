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
from datetime import date, timedelta
import glob
import scipy
import os
import netCDF4 as nc4
import time
import csv

#######################################
# iterate through time periods
#######################################
tf=1
for o in range(tf):
	startdate = date(2016,10,3)
# 	date1=date(starty[o],startm[o],1)
	date1=date(2017,1,1)
	delta = timedelta(90)
	offset = str(date1 - delta) 

	#######################################
	# read in the names of files to format
	# read in the names of .nc output
	#######################################
	with open ("LIST_STATION.csv") as myfile:
		datafile = myfile.read().split('\n')
	
	with open ("LIST_STATE.csv") as myfile:
		outputfile = myfile.read().split('\n')
	
	#######################################
	# define some time variables
	# Used for searching snoTEL data
	#######################################
	telstart='2016-10-03'
	telend='2017-01-01'

	#######################################
	# IMOPORT .csv files and find time .
	#######################################
	tt=len(datafile)
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
		
		#Some data formatting to remove any NaNs and make sure TMAX is > TMIN
		t_min = data['temperature_min'].values[start:finish]
		NaNs = np.isnan(t_min)
		length_NaNs=len(NaNs)
		for i in range(length_NaNs):
			nan = NaNs[i]
			if nan:
				t_min[i]=t_min[i-1]	
		
		t_max = data['temperature_max'].values[start:finish]
		NaNs = np.isnan(t_max)
		length_NaNs=len(NaNs)
		for i in range(length_NaNs):
			nan = NaNs[i]
			if nan:
				t_max[i]=t_max[i-1]
		
		prec = data['precipitation'].values[start:finish]
		NaNs = np.isnan(prec)
		length_NaNs=len(NaNs)
		for i in range(length_NaNs):
			nan = NaNs[i]
			if nan:
				prec[i]=prec[i-1]

		#######################################
		# netCDF creation
		#######################################
		outfile=outputfile[z]
		
		ncid = nc4.Dataset(outfile, "w", format="NETCDF4")

		dimid_T = ncid.createDimension('time',90)
		dimid_lat = ncid.createDimension('lat',1)
		dimid_lon = ncid.createDimension('lon',1)

		####################################### Variable: Time
		time_varid = ncid.createVariable('time','i8',('time',))
		length_time=90
		time2 = [i for i in range(length_time)]

		# Attributes
		time_varid.units         = 'days since 2016-10-03 00:00:00'
		time_varid.calendar      = 'proleptic_gregorian'

		# Write data
		time_varid[:] = time2

		####################################### Variables: Latitude & Longitude
		lon_varid = ncid.createVariable('lon','d',('lon',))
		lat_varid = ncid.createVariable('lat','d',('lat',))


		# Attributes
		lat_varid._fillvalue     = 'nan'
		lon_varid._fillvalue     = 'nan'
		lat_varid.standard_name  = 'latitude'
		lon_varid.standard_name  = 'longitude'
		lat_varid.units          = 'degrees_north'
		lon_varid.units          = 'degrees_east'
		lat_varid.axis           = 'Y'
		lon_varid.axis           = 'X'

		# Write data
		lat_varid[:] = lat
		lon_varid[:] = lon

		####################################### Variable: Precipitation
		prec_varid = ncid.createVariable('prec','d',('time','lat','lon',))

		# Attributes
		prec_varid.names         = '_fillvalue'
		prec_varid.value        = 'nan'

		# Write data
		prec_varid[:] = prec

		###################################### Variable: Temp minimum
		t_min_varid = ncid.createVariable('t_min','d',('time','lat','lon',))

		# Attributes
		t_min_varid.names         = '_fillvalue'
		t_min_varid.value        = 'nan'

		# Write data
		t_min_varid[:] = t_min

		####################################### Variable: Temp maximum
		t_max_varid = ncid.createVariable('t_max','d',('time','lat','lon',))

		# Attributes
		t_max_varid.names         = '_fillvalue'
		t_max_varid.value        = 'nan'

		# Write data
		t_max_varid[:] = t_max

		####################################### Header 
		ncid.License     = 'The file was created by C.Hart, https://github.com/ChristianHart2019'
		ncid.history     = 'Created ' + time.ctime(time.time())

		#####################
		ncid.close()


