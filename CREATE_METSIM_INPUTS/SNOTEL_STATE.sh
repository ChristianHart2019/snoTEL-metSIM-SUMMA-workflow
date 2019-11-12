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
# define some time variables
# Used for searching snoTEL data
#######################################
with open ("LIST_START.csv") as myfile:
    startfile = myfile.read().split('\n')	#Read list containing the end of the state data time-series, which is = to start of the simulation

timedata = pd.read_csv('LIST_TIME.csv', sep=',') # read in general information about the start and end of the simulation period

month = timedata['month'].values[:]
startm = timedata['startm'].values[:]
starty = timedata['starty'].values[:]
endd = timedata['endd'].values[:] 
endm = timedata['endm'].values[:] 
#######################################
# iterate through months
#######################################
tf=13
for o in range(tf):
	startdate = date(1900,1,1)
	date1=date(starty[o],startm[o],1)
	delta = timedelta(90)
	offset = str(date1 - delta) 
	
	telstart=offset	
	telend=startfile[o]

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
	telstart='2017-05-03'
	telend='2017-08-01'

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

		dimid_T = ncid.createDimension('time',90)
		dimid_lat = ncid.createDimension('lat',1)
		dimid_lon = ncid.createDimension('lon',1)

		####################################### Variable: Time
		time_varid = ncid.createVariable('time','i8',('time',))
		length_time=90
		time2 = [i for i in range(length_time)]

		# Attributes
		time_varid.units         = 'days since 1949-10-03 00:00:00'
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
	
		lat1 = [30.0312500000000, 30.0937500000000, 30.1562500000000, 30.2187500000000, 30.2812500000000, 30.3437500000000, 30.4062500000000, 30.4687500000000, 30.5312500000000]	
		lon1 = [-100.031250000000, -99.9687500000000, -99.9062500000000, -99.8437500000000, -99.7812500000000, -99.7187500000000, -99.6562500000000, -99.5937500000000, -99.5312500000000]
	
		length_lon=len(lon1)
		for i in range(length_lon): 
			lat1[i] = lat-(0.001*i)
			lon1[i] = lon-(0.001*i)
	
		# Write data
		lat_varid[:] = lat
		lon_varid[:] = lon

		####################################### Variable: Precipitation
		prec_varid = ncid.createVariable('prec','d',('time','lat','lon',))

		# Attributes
		prec_varid.names         = '_fillvalue'
		prec_varid.value        = 'nan'

		# Write data
		prec_length = len(prec)
		lat_length = len(lat1)
		lon_length = len(lon1)
		
		PREC = np.zeros((prec_length,lat_length,lon_length))	
		
		for i in range(prec_length):
			for j in range(lat_length):
				for k in range(lon_length):	
					PREC[i][j][k]=prec[i]

	# 	prec_varid[:] = PREC
		prec_varid[:] = prec

		###################################### Variable: Temp minimum
		t_min_varid = ncid.createVariable('t_min','d',('time','lat','lon',))

		# Attributes
		t_min_varid.names         = '_fillvalue'
		t_min_varid.value        = 'nan'

		# Write data
		t_min_length = len(t_min)
	
		T_MIN = np.zeros((t_min_length,lat_length,lon_length))
	
		for i in range(t_min_length):
			for j in range(lat_length):
				for k in range(lon_length):	
					T_MIN[i][j][k]=t_min[i]
	
	# 	t_min_varid[:] = T_MIN
		t_min_varid[:] = t_min

		####################################### Variable: Temp maximum
		t_max_varid = ncid.createVariable('t_max','d',('time','lat','lon',))

		# Attributes
		t_max_varid.names         = '_fillvalue'
		t_max_varid.value        = 'nan'

		# Write data
		t_max_length = len(t_max)
	
		T_MAX = np.zeros((t_max_length,lat_length,lon_length))
	# 	
		for i in range(t_max_length):
			for j in range(lat_length):
				for k in range(lon_length):	
					T_MAX[i][j][k]=t_max[i]
	
	# 	t_max_varid[:] = T_MAX
		t_max_varid[:] = t_min

		####################################### Variable: Snow water equivelent
	# 	swe_varid = ncid.createVariable('swe','d',('time',))
	# 
	# 	# Attributes
	# 	swe_varid.names         = '_fillvalue'
	# 	swe_varid.value        = 'nan'
	# 
	# 	# Write data
	# 	swe_varid[:] = swe

		####################################### Header 
		ncid.License     = 'The file was created by C.Hart, https://github.com/ChristianHart2019'
		ncid.history     = 'Created ' + time.ctime(time.time())

		#####################
		ncid.close()


