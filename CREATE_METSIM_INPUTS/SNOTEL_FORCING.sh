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
	
	with open ("LIST_FORCING.csv") as myfile:
		outputfile = myfile.read().split('\n')
	
	#######################################
	# define some time variables
	# Used for searching snoTEL data
	# must be 1 day longer than the config
	#######################################
# 	telstart='2017-08-01'
# 	telend='2017-09-01'

	telstart=str(startfile[o])
	telend=str(endfile[o])
	
# 	print(telstart)
# 	print(telend)	
	# telstart='1950-01-01'
	# telend='1950-01-31'

	comj = (endd[o]==31)  			#Calculate the number of days in each month and
	if comj:						#how they correspond to the configuration file
		endday=1
		endmonth=2
	coma = (endd[o]==30)
	if coma:
		endday=31
		endmonth=1
	comb = (endd[o]==28)
	if comb:
		endday=29
		endmonth=1
	
	#######################################
	# The date range used to declare this variable has to be 1 day
	# longer than the range used in the CONFIG file, but it is not the actual 
	# range of the data. For MetSim to work this range will have to be from a time period
	# originally declared in the config file starting in 1950/1/1
	########################################
	startdate = date(1900,1,1)
	date1 = date(1950,1,1)
	date2 = date(1950,endmonth,endday)
# 	date1 = date(1950,8,1)
# 	date2 = date(1950,10,1)
	time1 = (date1-startdate).days		
	time2 = (date2-startdate).days		# Calculate days since 1900 for the start and end of time-series
	td=time2-time1						# Calculate length of output data
	int_time1 = int(time1)
	int_time2 = int(time2)

	time3 = np.arange(int_time1, int_time2)
	
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
			
		print(start)
		print(finish)
				
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
		Tmin = data['temperature_min'].values[start:finish]
		Tmax = data['temperature_max'].values[start:finish]
		Prec = data['precipitation'].values[start:finish]

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
		dimid_T = ncid.createDimension('time',td)

		####################################### Variables: Latitude & Longitude
		lon_varid = ncid.createVariable('lon','d',('lon',))
		lat_varid = ncid.createVariable('lat','d',('lat',))

		# Attributes
		lat_varid.standard_name  = 'latitude'
		lon_varid.standard_name  = 'longitude'
		lat_varid.long_name      = 'latitude'
		lon_varid.long_name      = 'longitude'
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
	
		####################################### Variable: Time
		### The date range used to declare this variable has to be 1 day
		### longer than the range used in the CONFIG file, but it is not the actual 
		### range of the data. For MetSim to work this range will have to be from a time period
		### originally declared in the config file starting in 1950/1/1
		########################################
		time_varid = ncid.createVariable('time','d',('time',))

		# Attributes
		time_varid.standard_name  = 'time'
		time_varid.longname       = 'Time axis'
		time_varid.units          = 'days since 1900-01-01 00:00:00'
		time_varid.calendar       = 'standard'

		# Write data
		time_varid[:] = time3

		####################################### Variable: Precipitation
		Prec_varid = ncid.createVariable('Prec','f',('time','lat','lon',))
	# 	Prec_varid = ncid.createVariable('Prec','f',('time',))

		# Attributes
		Prec_varid.units            = 'mm'
		Prec_varid.FillValue        = '1.00e+20'
		Prec_varid.missingvalue     = '1.00e+20'
		Prec_varid.long_name        = 'precipitation'

		# Write data
		Prec_length = len(Prec)
		lat_length = len(lat1)
		lon_length = len(lon1)
		
		PREC = np.zeros((Prec_length,lat_length,lon_length))	
		
		for i in range(Prec_length):
			for j in range(lat_length):
				for k in range(lon_length):	
					PREC[i][j][k]=Prec[i]+0.001		# had to add a small number here to get the correct decimal place

	# 	Prec_varid[:] = PREC
		Prec_varid[:] = Prec
		
		####################################### Variable: Temp maximum
		Tmax_varid = ncid.createVariable('Tmax','f',('time','lat','lon',))
	# 	Tmax_varid = ncid.createVariable('Tmax','f',('time',))

		# Attributes
		Tmax_varid.units            = 'C'
		Tmax_varid.FillValue        = '1.00e+20'
		Tmax_varid.missingvalue    = '1.00e+20'
		Tmax_varid.long_name        = 'Daily maximum temperature'

		# Write data
		Tmax_length = len(Tmax)
		lat_length = len(lat1)
		lon_length = len(lon1)
		
		TMAX = np.zeros((Tmax_length,lat_length,lon_length))	
		
		for i in range(Tmax_length):
			for j in range(lat_length):
				for k in range(lon_length):	
					TMAX[i][j][k]=Tmax[i]

	# 	Tmax_varid[:] = TMAX
		Tmax_varid[:] = Tmax
		
		###################################### Variable: Temp minimum
		Tmin_varid = ncid.createVariable('Tmin','f',('time','lat','lon',))
	# 	Tmin_varid = ncid.createVariable('Tmin','f',('time',))

		# Attributes
		Tmin_varid.units            = 'C'
		Tmin_varid.FillValue        = '1.00e+20'
		Tmin_varid.missingvalue    = '1.00e+20'
		Tmin_varid.long_name        = 'Daily minimum temperature'

		# Write data
		Tmin_length = len(Tmin)
		lat_length = len(lat1)
		lon_length = len(lon1)
		
		TMIN = np.zeros((Tmin_length,lat_length,lon_length))	
		
		for i in range(Tmin_length):
			for j in range(lat_length):
				for k in range(lon_length):	
					TMIN[i][j][k]=Tmin[i]

	# 	Tmin_varid[:] = TMIN
		Tmin_varid[:] = Tmin

		####################################### Header 
		ncid.License     = 'The file was created by C.Hart, https://github.com/ChristianHart2019'
		ncid.history     = 'Created ' + time.ctime(time.time())

		#####################
		ncid.close()


