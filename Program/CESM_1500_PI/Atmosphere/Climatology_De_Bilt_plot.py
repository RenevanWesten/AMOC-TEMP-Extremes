#Program plots the climatology for De Bilt

from pylab import *
import numpy
import time
import glob, os
import math
import netCDF4 as netcdf
from scipy import stats
from cartopy import crs as ccrs, feature as cfeature
import matplotlib.ticker as mticker
from scipy.stats import genextreme

#Making pathway to folder with all data
directory_1500_PI = '../../../Data/CESM_1500_PI/'
directory_obs	  = '../../../Data/Observations/'

def ReturnValue(prob, shape, loc = 0.0, scale = 1.0):
    """Return the return value at a given probability"""
	
    return loc - (scale / shape)* (1.0 - (-np.log(1 - prob))**(-shape))

def ReturnTime(value, shape, loc = 0.0, scale = 1.0):
    """Returns the return time of a given event"""

    prob	= 1.0 - np.exp(-(1.0 + shape * ( (value - loc) / scale))**(-1.0 / shape))

    return 1.0 / prob

def wind_uv_to_dir(U,V):
    """Calculates the wind direction from the u and v component of wind.
    Takes into account the wind direction coordinates is different than the 
    trig unit circle coordinate. If the wind directin is 360 then returns zero (by %360)
    Inputs:
        U = west/east direction (wind from the west is positive, from the east is negative)
        V = south/north direction (wind from the south is positive, from the north is negative)"""
	
    return (270-np.rad2deg(np.arctan2(V,U)))%360
	
#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

fh = netcdf.Dataset(directory_1500_PI+'Atmosphere/Climatology_De_Bilt.nc', 'r')

#Writing data to correct variable
time_all	   = fh.variables['time'][:] 			
day_cesm	   = fh.variables['day'][:] 			
temp_all	   = fh.variables['TEMP'][:] 			
temp_all_min   = fh.variables['TEMP_min'][:] 			
temp_all_max   = fh.variables['TEMP_max'][:] 			
pres_all	   = fh.variables['PSL'][:] 	
u_vel_all	   = fh.variables['U_PSL'][:]
v_vel_all	   = fh.variables['V_PSL'][:]

fh.close()

#-----------------------------------------------------------------------------------------

wind_dir 		        = wind_uv_to_dir(u_vel_all, v_vel_all)
wind_dir[wind_dir>180]	= wind_dir[wind_dir>180] - 360.0
index_min		        = np.argmin(temp_all_min, axis = 1)

wind_dir_cold_extreme	= np.zeros(len(time_all))
wind_speed_cold_extreme	= np.zeros(len(time_all))
pres_cold_extreme	    = np.zeros(len(time_all))
wind_dir_warm_extreme	= np.zeros(len(time_all))
wind_speed_warm_extreme	= np.zeros(len(time_all))
pres_warm_extreme	    = np.zeros(len(time_all))

temp_winter		       = np.zeros((len(time_all), 31+31+28))
wind_dir_winter		   = np.zeros((len(time_all), 31+31+28))
temp_summer		       = np.zeros((len(time_all), 30+31+31))
wind_dir_summer		   = np.zeros((len(time_all), 30+31+31))

for time_i in range(len(time_all)):
	#Determine the wind direction for the various seasons
	wind_dir_cold_extreme[time_i]	= wind_dir[time_i, np.argmin(temp_all[time_i])]
	wind_speed_cold_extreme[time_i]	= np.sqrt(u_vel_all[time_i, np.argmin(temp_all[time_i])]**2 + u_vel_all[time_i, np.argmin(temp_all[time_i])]**2)
	pres_cold_extreme[time_i]	    = pres_all[time_i, np.argmin(temp_all[time_i])]
	wind_dir_warm_extreme[time_i]	= wind_dir[time_i, np.argmax(temp_all[time_i])]
	wind_speed_warm_extreme[time_i]	= np.sqrt(u_vel_all[time_i, np.argmax(temp_all[time_i])]**2 + u_vel_all[time_i, np.argmax(temp_all[time_i])]**2)
	pres_warm_extreme[time_i]	    = pres_all[time_i, np.argmax(temp_all[time_i])]
	
	temp_winter[time_i, :31+28]	    = temp_all[time_i, :31+28]	#January+February, year n
	temp_winter[time_i, 31+28:]	    = temp_all[time_i, -31:]	#December, year n
	wind_dir_winter[time_i, :31+28]	= wind_dir[time_i, :31+28]	#January+February, year n
	wind_dir_winter[time_i, 31+28:]	= wind_dir[time_i, -31:]	#December, year n	
	
	temp_summer[time_i]		        = temp_all[time_i, 31+28+31+30+31:31+28+31+30+31+30+31+31]	#June+July+August, year n
	wind_dir_summer[time_i]		    = wind_dir[time_i, 31+28+31+30+31:31+28+31+30+31+30+31+31]	#June+July+August, year n	

temp_winter		= temp_winter.ravel()
wind_dir_winter	= wind_dir_winter.ravel()
temp_summer		= temp_summer.ravel()
wind_dir_summer	= wind_dir_summer.ravel()

#-----------------------------------------------------------------------------------------

#Make bins with the different percentiles and for the different (8) wind direction
temp_winter_hist	= np.zeros((5, 8))

temp_winter_1		= np.percentile(temp_winter, 1)
temp_winter_5		= np.percentile(temp_winter, 5)
temp_winter_10		= np.percentile(temp_winter, 10)
temp_winter_25		= np.percentile(temp_winter, 25)

for day_i in range(len(temp_winter)):
    #Get the wind direction
    wind_dir_day	= wind_dir_winter[day_i]
	
    if wind_dir_day < -22.5: 
        #In its form to start with north as the first index
        wind_dir_day	= wind_dir_day + 360.0
		
    wind_dir_index	= int((wind_dir_day + 22.5) / 45)
	
    if temp_winter[day_i] <= temp_winter_1:
        #Lower than 1%
        temp_winter_hist[0, wind_dir_index] += 1

    elif temp_winter[day_i] > temp_winter_1 and temp_winter[day_i] <= temp_winter_5:
        #1% to 5%
        temp_winter_hist[1, wind_dir_index] += 1

    elif temp_winter[day_i] > temp_winter_5 and temp_winter[day_i] <= temp_winter_10:
        #5% to 10%
        temp_winter_hist[2, wind_dir_index] += 1
		
    elif temp_winter[day_i] > temp_winter_10 and temp_winter[day_i] <= temp_winter_25:
        #10% to 25%
        temp_winter_hist[3, wind_dir_index] += 1

    elif temp_winter[day_i] > temp_winter_25:
        #Above 25%
        temp_winter_hist[4, wind_dir_index] += 1
		
print('Days with a NW component for 10% coldest days:', np.sum(temp_winter_hist[:3, -1]))

for temp_i in range(len(temp_winter_hist)):
    #Now norm for each temperature bin
    temp_winter_hist[temp_i]	= temp_winter_hist[temp_i] / np.sum(temp_winter_hist[temp_i])

#-----------------------------------------------------------------------------------------

#Make bins with the different percentiles and for the different (8) wind direction
temp_summer_hist	= np.zeros((5, 8))

temp_summer_75		= np.percentile(temp_summer, 75)
temp_summer_90		= np.percentile(temp_summer, 90)
temp_summer_95		= np.percentile(temp_summer, 95)
temp_summer_99		= np.percentile(temp_summer, 99)

for day_i in range(len(temp_summer)):
    #Get the wind direction
    wind_dir_day	= wind_dir_summer[day_i]
	
    if wind_dir_day < -22.5: 
        #In its form to start with north as the first index
        wind_dir_day	= wind_dir_day + 360.0
		
    wind_dir_index	= int((wind_dir_day + 22.5) / 45)
	
    if temp_summer[day_i] >= temp_summer_99:
        #Higher than 99%
        temp_summer_hist[4, wind_dir_index] += 1

    elif temp_summer[day_i] >= temp_summer_95 and temp_summer[day_i] < temp_summer_99:
        #95% to 99%
        temp_summer_hist[3, wind_dir_index] += 1

    elif temp_summer[day_i] >= temp_summer_90 and temp_summer[day_i] < temp_summer_95:
        #90% to 95%
        temp_summer_hist[2, wind_dir_index] += 1
		
    elif temp_summer[day_i] >= temp_summer_75 and temp_summer[day_i] < temp_summer_90:
        #75% to 95%
        temp_summer_hist[1, wind_dir_index] += 1

    elif temp_summer[day_i] < temp_summer_75:
        #Below 75%
        temp_summer_hist[0, wind_dir_index] += 1


for temp_i in range(len(temp_summer_hist)):
    #Now norm for each temperature bin
    temp_summer_hist[temp_i]	= temp_summer_hist[temp_i] / np.sum(temp_summer_hist[temp_i])

#Convert back to 0 - 360E
wind_dir_warm_extreme[wind_dir_warm_extreme<0]	= wind_dir_warm_extreme[wind_dir_warm_extreme<0] + 360.0

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.plot(np.min(temp_all, axis = 1), wind_dir_cold_extreme, 'o', color = 'k')
ax.set_xlim(-41, 11)
ax.set_ylim(-190, 190)

ax.set_xlabel('Daily-averaged temperature ($^{\circ}$C)')
ax.set_ylabel('Wind direction')
ax.set_yticks(np.arange(-180, 180.1, 45))
ax.set_yticklabels(['S', 'SW', 'W', 'NW', 'N', 'NE', 'E', 'SE', 'S'])
ax.grid()
ax.set_title('c) Cold extreme and wind direction, AMOC on ($F_H$ = 0.45 Sv)')

#-----------------------------------------------------------------------------------------

ax2 	= fig.add_axes([0.20, 0.43, 0.28, 0.40])

wind_label = ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW']
wind_color = ['darkviolet', 'deepskyblue', 'royalblue', 'pink', 'firebrick', 'sandybrown', 'khaki', 'limegreen']

bottom	= np.zeros(len(temp_winter_hist))

for wind_i in range(8):	

	ax2.bar(np.arange(0, 4.1, 1), -temp_winter_hist[:, wind_i], 0.8, label = wind_label[wind_i], color = wind_color[wind_i], bottom = 1.0-bottom)
	bottom	+= temp_winter_hist[:, wind_i]

ax2.legend(loc = (1.02, -0.50), framealpha = 1)

ax2.set_xlim(-0.5, 4.5)
ax2.set_ylim(0, 1)

ax2.set_xticks(np.arange(0, 5))
ax2.set_xticklabels(['$\leq 1\%$ ($'+str(round(temp_winter_1,1))+'^{\circ}$C)', '$1\%$ to $5\%$ ($'+str(round(temp_winter_5,1))+'^{\circ}$C)', '$5\%$ to $10\%$ ($'+str(round(temp_winter_10,1))+'^{\circ}$C)', '$10\%$ to $25\%$ ($'+str(round(temp_winter_25,1))+'^{\circ}$C)', '$> 25\%$'], rotation=-90, fontsize = 9)

ax2.set_title('Winter (DJF)')

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.plot(np.max(temp_all, axis = 1), wind_dir_warm_extreme, 'o', color = 'k')
ax.set_xlim(14, 51)
ax.set_ylim(-10, 370)

ax.set_xlabel('Daily-averaged temperature ($^{\circ}$C)')
ax.set_ylabel('Wind direction')
ax.set_yticks(np.arange(0, 360.1, 45))
ax.set_yticklabels(['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW', 'N'])
ax.grid()
ax.set_title('c) Warm extreme and direction, AMOC on ($F_H$ = 0.45 Sv)')

#-----------------------------------------------------------------------------------------

ax2 	= fig.add_axes([0.58, 0.43, 0.28, 0.40])

wind_label = ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW']
wind_color = ['darkviolet', 'deepskyblue', 'royalblue', 'pink', 'firebrick', 'sandybrown', 'khaki', 'limegreen']

bottom	= np.zeros(len(temp_summer_hist))

for wind_i in range(8):	

	ax2.bar(np.arange(0, 4.1, 1), -temp_summer_hist[:, wind_i], 0.8, label = wind_label[wind_i], color = wind_color[wind_i], bottom = 1.0-bottom)
	bottom	+= temp_summer_hist[:, wind_i]

ax2.legend(loc = (1.02, 0.05), framealpha = 1)

ax2.set_xlim(-0.5, 4.5)
ax2.set_ylim(0, 1)

ax2.set_xticks(np.arange(0, 5))
ax2.set_xticklabels(['$< 75\%$ ($'+str(round(temp_summer_75,1))+'^{\circ}$C)', '$75\%$ to $90\%$ ($'+str(round(temp_summer_90,1))+'^{\circ}$C)', '$90\%$ to $95\%$ ($'+str(round(temp_summer_95,1))+'^{\circ}$C)', '$95\%$ to $99\%$ ($'+str(round(temp_summer_99,1))+'^{\circ}$C)', '$\geq 99\%$'], rotation=-90, fontsize = 9)

ax2.set_title('Summer (JJA)')

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

#Get the GEV fits
temp_min_extreme_1			= -np.min(temp_all_min, axis = 1)
temp_min_extreme_2			= -np.min(temp_all, axis = 1)
temp_min_extreme_3			= -np.min(temp_all_max, axis = 1)

shape_min_1, loc_min_1, scale_min_1	= genextreme.fit(temp_min_extreme_1, loc = np.mean(temp_min_extreme_1), scale = np.std(temp_min_extreme_1))
shape_min_2, loc_min_2, scale_min_2	= genextreme.fit(temp_min_extreme_2, loc = np.mean(temp_min_extreme_2), scale = np.std(temp_min_extreme_2))
shape_min_3, loc_min_3, scale_min_3	= genextreme.fit(temp_min_extreme_3, loc = np.mean(temp_min_extreme_3), scale = np.std(temp_min_extreme_3))
shape_min_1, shape_min_2, shape_min_3	= -shape_min_1, -shape_min_2, -shape_min_3

fig, ax = subplots()

prob		= np.linspace(1.0 / 100.0, 0.9999, 100000)					#Fit
z_min_1    	= ReturnValue(prob, shape_min_1, -loc_min_1, -scale_min_1) #Fit curve
z_min_2    	= ReturnValue(prob, shape_min_2, -loc_min_2, -scale_min_2) #Fit curve
z_min_3    	= ReturnValue(prob, shape_min_3, -loc_min_3, -scale_min_3) #Fit curve

freq 		= np.linspace(1.0 / len(temp_min_extreme_1), 1, len(temp_min_extreme_1))	#Observed frequency

ax.plot(1.0/freq , sorted(-temp_min_extreme_1), 'o', alpha = 0.5, color = 'royalblue')
ax.plot(1.0/freq , sorted(-temp_min_extreme_2), 'o', alpha = 0.5, color = 'k')
ax.plot(1.0/freq , sorted(-temp_min_extreme_3), 'o', alpha = 0.5, color = 'firebrick')

graph_3    = ax.plot(1.0 / prob, z_min_3, '-', color = 'firebrick', linewidth = 2.0, label = 'Daily maxima')
graph_2    = ax.plot(1.0 / prob, z_min_2, '-', color = 'k', linewidth = 2.0, label = 'Daily average')
graph_1    = ax.plot(1.0 / prob, z_min_1, '-', color = 'royalblue', linewidth = 2.0, label = 'Daily minima')

ax.set_xlim(1, 100)
ax.set_ylim(-41, 11)
ax.grid()
ax.set_xlabel('Return time (years)')
ax.set_ylabel('Temperature ($^{\circ}$C)')
ax.set_xscale('log')
ax.set_xticks([1, 3, 5, 10, 30, 50, 100])
ax.set_xticklabels([1, 3, 5, 10, 30, 50, 100])

ax.legend(loc = 'lower left', framealpha = 1)

ax.set_title('d) Cold extreme return times, AMOC on ($F_H$ = 0.45 Sv)')

print()
print('1:10-year cold extreme for daily minimum:', ReturnValue(1/10, shape_min_1, -loc_min_1, -scale_min_1))
print('1:10-year cold extreme for daily average:', ReturnValue(1/10, shape_min_2, -loc_min_2, -scale_min_2)) 
print('1:10-year cold extreme for daily maximum:', ReturnValue(1/10, shape_min_3, -loc_min_3, -scale_min_3)) 

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

temp_max_extreme_1			= np.max(temp_all_min, axis = 1)
temp_max_extreme_2			= np.max(temp_all, axis = 1)
temp_max_extreme_3			= np.max(temp_all_max, axis = 1)

shape_max_1, loc_max_1, scale_max_1	= genextreme.fit(temp_max_extreme_1, loc = np.mean(temp_max_extreme_1), scale = np.std(temp_max_extreme_1))
shape_max_2, loc_max_2, scale_max_2	= genextreme.fit(temp_max_extreme_2, loc = np.mean(temp_max_extreme_2), scale = np.std(temp_max_extreme_2))
shape_max_3, loc_max_3, scale_max_3	= genextreme.fit(temp_max_extreme_3, loc = np.mean(temp_max_extreme_3), scale = np.std(temp_max_extreme_3))
shape_max_1, shape_max_2, shape_max_3	= -shape_max_1, -shape_max_2, -shape_max_3


fig, ax = subplots()

prob		= np.linspace(1.0 / 100.0, 0.9999, 100000)					#Fit
z_max_1    	= ReturnValue(prob, shape_max_1, loc_max_1, scale_max_1) 	#Fit curve
z_max_2    	= ReturnValue(prob, shape_max_2, loc_max_2, scale_max_2) 	#Fit curve
z_max_3    	= ReturnValue(prob, shape_max_3, loc_max_3, scale_max_3) 	#Fit curve

freq 		= np.linspace(1.0 / len(temp_max_extreme_1), 1, len(temp_max_extreme_1))	#Observed frequency

ax.plot(1.0/freq , sorted(temp_max_extreme_1)[::-1], 'o', alpha = 0.5, color = 'royalblue')
ax.plot(1.0/freq , sorted(temp_max_extreme_2)[::-1], 'o', alpha = 0.5, color = 'k')
ax.plot(1.0/freq , sorted(temp_max_extreme_3)[::-1], 'o', alpha = 0.5, color = 'firebrick')

graph_3    = ax.plot(1.0 / prob, z_max_3, '-', color = 'firebrick', linewidth = 2.0, label = 'Daily maxima')
graph_2    = ax.plot(1.0 / prob, z_max_2, '-', color = 'k', linewidth = 2.0, label = 'Daily average')
graph_1    = ax.plot(1.0 / prob, z_max_1, '-', color = 'royalblue', linewidth = 2.0, label = 'Daily minima')

ax.set_xlim(1, 100)
ax.set_ylim(14, 51)
ax.grid()
ax.set_xlabel('Return time (years)')
ax.set_ylabel('Temperature ($^{\circ}$C)')
ax.set_xscale('log')
ax.set_xticks([1, 3, 5, 10, 30, 50, 100])
ax.set_xticklabels([1, 3, 5, 10, 30, 50, 100])

ax.legend(loc = 'lower right', framealpha = 1)

ax.set_title('d) Warm extreme return times, AMOC on ($F_H$ = 0.45 Sv)')

print()
print('1:10-year warm extreme for daily minimum:', ReturnValue(1/10, shape_max_1, loc_max_1, scale_max_1))
print('1:10-year warm extreme for daily average:', ReturnValue(1/10, shape_max_2, loc_max_2, scale_max_2)) 
print('1:10-year warm extreme for daily maximum:', ReturnValue(1/10, shape_max_3, loc_max_3, scale_max_3))

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

temp_obs		= ma.masked_all((100, 365))
temp_obs_min	= ma.masked_all((100, 365))
temp_obs_max	= ma.masked_all((100, 365))

#Read in the data from De Bilt
fh      = open(directory_obs+'DeBilt.txt', 'r')
lines   = fh.readlines()
fh.close()

year_start	= 1901
day_counter	= 0

for line_i in range(53, len(lines)):
    #Read in the data
    line 	= lines[line_i].split(',')
    date	= str(line[1])
    year	= int(date[0:4])
    month	= int(date[4:6])
    day	= int(date[6:8])

    if year > 2000:
        #Only consider 1901 - 2000
        break

    if month == 2 and day == 29:
        #No leap year
        continue

    if year != year_start:
        #Happy new year! Update counters
        year_start 	+= 1
        day_counter	= 0

    #Get the daily temperature statistics (in degC)
    temp		= int(line[11]) / 10.0
    temp_min	= int(line[12]) / 10.0
    temp_max	= int(line[14]) / 10.0

    temp_obs[year-1901, day_counter]		= temp
    temp_obs_min[year-1901, day_counter]	= temp_min
    temp_obs_max[year-1901, day_counter]	= temp_max

    day_counter += 1

#-----------------------------------------------------------------------------------------

ice_days_ref 		= round(len(np.where(temp_obs_max < 0.0)[0]) / len(temp_obs), 1)
frost_days_ref 		= round(len(np.where(temp_obs_min < 0.0)[0]) / len(temp_obs), 1)
warm_days_ref 		= round(len(np.where(temp_obs_max >= 20.0)[0]) / len(temp_obs), 1)
summer_days_ref 	= round(len(np.where(temp_obs_max >= 25.0)[0]) / len(temp_obs), 1)
tropical_days_ref 	= round(len(np.where(temp_obs_max >= 30.0)[0]) / len(temp_obs), 1)
tropical_nights_ref	= round(len(np.where(temp_obs_min >= 20.0)[0]) / len(temp_obs), 1)

ice_days		    = round(len(np.where(temp_all_max < 0.0)[0]) / len(temp_all), 1)
frost_days 		    = round(len(np.where(temp_all_min < 0.0)[0]) / len(temp_all), 1)
warm_days 		    = round(len(np.where(temp_all_max >= 20.0)[0]) / len(temp_all), 1)
summer_days 		= round(len(np.where(temp_all_max >= 25.0)[0]) / len(temp_all), 1)
tropical_days		= round(len(np.where(temp_all_max >= 30.0)[0]) / len(temp_all), 1)
tropical_nights		= round(len(np.where(temp_all_min >= 20.0)[0]) / len(temp_all), 1)


#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.fill_between(day_cesm, np.min(temp_obs, axis = 0), np.max(temp_obs, axis = 0), alpha=0.25, edgecolor='mediumorchid', facecolor='mediumorchid')
ax.plot(day_cesm, np.mean(temp_obs, axis = 0), '-', color = 'mediumorchid', label = 'Observations')
ax.plot(day_cesm, np.max(temp_all, axis = 0), '-', color = 'firebrick', label = 'Highest')
ax.plot(day_cesm, np.mean(temp_all, axis = 0), '-', color = 'k', label = 'Mean')
ax.plot(day_cesm, np.min(temp_all, axis = 0), '-', color = 'royalblue', label = 'Lowest')

ax.set_xlim(1, 365)
ax.set_ylim(-45, 45)
ax.set_ylabel('Temperature ($^{\circ}$C)')
ax.grid()

ax.set_xticks([16, 46, 76, 106, 136, 167, 197, 228, 259, 289, 320, 350], minor=False)
ax.set_xticks([1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335], minor=True)
ax.set_xticklabels(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
ax.xaxis.grid(False, which='major')
ax.xaxis.grid(True, which='minor')

ax.set_title('c) Daily average, AMOC on ($F_H$ = 0.45 Sv)')
legend_1	= ax.legend(loc=(0.265, 0.02), ncol=1, framealpha = 1.0, numpoints = 1)

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.fill_between(day_cesm, np.min(temp_obs_min, axis = 0), np.max(temp_obs_min, axis = 0), alpha=0.25, edgecolor='mediumorchid', facecolor='mediumorchid')
ax.plot(day_cesm, np.mean(temp_obs_min, axis = 0), '-', color = 'mediumorchid', label = 'Observations')	
ax.plot(day_cesm, np.max(temp_all_min, axis = 0), '-', color = 'firebrick', label = 'Highest')
ax.plot(day_cesm, np.mean(temp_all_min, axis = 0), '-', color = 'k', label = 'Mean')
ax.plot(day_cesm, np.min(temp_all_min, axis = 0), '-', color = 'royalblue', label = 'Lowest')

ax.set_xlim(1, 365)
ax.set_ylim(-45, 45)
ax.set_ylabel('Temperature ($^{\circ}$C)')
ax.grid()

ax.set_xticks([16, 46, 76, 106, 136, 167, 197, 228, 259, 289, 320, 350], minor=False)
ax.set_xticks([1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335], minor=True)
ax.set_xticklabels(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
ax.xaxis.grid(False, which='major')
ax.xaxis.grid(True, which='minor')

ax.set_title('e) Daily minima, AMOC on ($F_H$ = 0.45 Sv)')
legend_1	= ax.legend(loc=(0.265, 0.02), ncol=1, framealpha = 1.0, numpoints = 1)

ax.text(215, -22.5, '$T^{\mathrm{min}} < 0^{\circ}$C  : '+str(frost_days)+' yr$^{-1}$', color = 'darkturquoise', fontsize = 11, horizontalalignment='left', verticalalignment='center', zorder = 10)
ax.text(215, -27.5, '$T^{\mathrm{min}} \geq 20^{\circ}$C: '+str(tropical_nights)+' yr$^{-1}$', color = 'deeppink', fontsize = 11, horizontalalignment='left', verticalalignment='center', zorder = 10)

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.fill_between(day_cesm, np.min(temp_obs_max, axis = 0), np.max(temp_obs_max, axis = 0), alpha=0.25, edgecolor='mediumorchid', facecolor='mediumorchid')
ax.plot(day_cesm, np.mean(temp_obs_max, axis = 0), '-', color = 'mediumorchid', label = 'Observations')	
ax.plot(day_cesm, np.max(temp_all_max, axis = 0), '-', color = 'firebrick', label = 'Highest')
ax.plot(day_cesm, np.mean(temp_all_max, axis = 0), '-', color = 'k', label = 'Mean')
ax.plot(day_cesm, np.min(temp_all_max, axis = 0), '-', color = 'royalblue', label = 'Lowest')

ax.set_xlim(1, 365)
ax.set_ylim(-45, 45)
ax.set_ylabel('Temperature ($^{\circ}$C)')
ax.grid()

ax.set_xticks([16, 46, 76, 106, 136, 167, 197, 228, 259, 289, 320, 350], minor=False)
ax.set_xticks([1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335], minor=True)
ax.set_xticklabels(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
ax.xaxis.grid(False, which='major')
ax.xaxis.grid(True, which='minor')

ax.set_title('f) Daily maxima, AMOC on ($F_H$ = 0.45 Sv)')
legend_1	= ax.legend(loc=(0.265, 0.02), ncol=1, framealpha = 1.0, numpoints = 1)

ax.text(215, -22.5, '$T^{\mathrm{max}} < 0^{\circ}$C  : '+str(ice_days)+' yr$^{-1}$', color = 'dodgerblue', fontsize = 11, horizontalalignment='left', verticalalignment='center', zorder = 10)
ax.text(215, -27.5, '$T^{\mathrm{max}} \geq 20^{\circ}$C: '+str(warm_days)+' yr$^{-1}$', color = 'yellowgreen', fontsize = 11, horizontalalignment='left', verticalalignment='center', zorder = 10)
ax.text(215, -32.5, '$T^{\mathrm{max}} \geq 25^{\circ}$C: '+str(summer_days)+' yr$^{-1}$', color = 'sandybrown', fontsize = 11, horizontalalignment='left', verticalalignment='center', zorder = 10)
ax.text(215, -37.5, '$T^{\mathrm{max}} \geq 30^{\circ}$C: '+str(tropical_days)+' yr$^{-1}$', color = 'indianred', fontsize = 11, horizontalalignment='left', verticalalignment='center', zorder = 10)

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.plot(day_cesm, np.max(temp_obs, axis = 0), '-', color = 'firebrick', label = 'Highest')
ax.plot(day_cesm, np.mean(temp_obs, axis = 0), '-', color = 'k', label = 'Mean')
ax.plot(day_cesm, np.min(temp_obs, axis = 0), '-', color = 'royalblue', label = 'Lowest')

ax.set_xlim(1, 365)
ax.set_ylim(-45, 45)
ax.set_ylabel('Temperature ($^{\circ}$C)')
ax.grid()

ax.set_xticks([16, 46, 76, 106, 136, 167, 197, 228, 259, 289, 320, 350], minor=False)
ax.set_xticks([1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335], minor=True)
ax.set_xticklabels(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
ax.xaxis.grid(False, which='major')
ax.xaxis.grid(True, which='minor')

ax.set_title('a) Daily average, Observations (De Bilt, 1901 - 2000)')
legend_1	= ax.legend(loc=(0.265, 0.02), ncol=1, framealpha = 1.0, numpoints = 1)

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.plot(day_cesm, np.max(temp_obs_min, axis = 0), '-', color = 'firebrick', label = 'Highest')
ax.plot(day_cesm, np.mean(temp_obs_min, axis = 0), '-', color = 'k', label = 'Mean')
ax.plot(day_cesm, np.min(temp_obs_min, axis = 0), '-', color = 'royalblue', label = 'Lowest')

ax.set_xlim(1, 365)
ax.set_ylim(-45, 45)
ax.set_ylabel('Temperature ($^{\circ}$C)')
ax.grid()

ax.set_xticks([16, 46, 76, 106, 136, 167, 197, 228, 259, 289, 320, 350], minor=False)
ax.set_xticks([1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335], minor=True)
ax.set_xticklabels(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
ax.xaxis.grid(False, which='major')
ax.xaxis.grid(True, which='minor')

ax.set_title('a) Daily minima, Observations (De Bilt, 1901 - 2000)')
legend_1	= ax.legend(loc=(0.265, 0.02), ncol=1, framealpha = 1.0, numpoints = 1)

ax.text(215, -22.5, '$T^{\mathrm{min}} < 0^{\circ}$C  : '+str(frost_days_ref)+' yr$^{-1}$', color = 'darkturquoise', fontsize = 11, horizontalalignment='left', verticalalignment='center', zorder = 10)
ax.text(215, -27.5, '$T^{\mathrm{min}} \geq 20^{\circ}$C: '+str(tropical_nights_ref)+' yr$^{-1}$', color = 'deeppink', fontsize = 11, horizontalalignment='left', verticalalignment='center', zorder = 10)

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.plot(day_cesm, np.max(temp_obs_max, axis = 0), '-', color = 'firebrick', label = 'Highest')
ax.plot(day_cesm, np.mean(temp_obs_max, axis = 0), '-', color = 'k', label = 'Mean')
ax.plot(day_cesm, np.min(temp_obs_max, axis = 0), '-', color = 'royalblue', label = 'Lowest')

ax.set_xlim(1, 365)
ax.set_ylim(-45, 45)
ax.set_ylabel('Temperature ($^{\circ}$C)')
ax.grid()

ax.set_xticks([16, 46, 76, 106, 136, 167, 197, 228, 259, 289, 320, 350], minor=False)
ax.set_xticks([1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335], minor=True)
ax.set_xticklabels(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
ax.xaxis.grid(False, which='major')
ax.xaxis.grid(True, which='minor')

ax.set_title('b) Daily maxima, Observations (De Bilt, 1901 - 2000)')
legend_1	= ax.legend(loc=(0.265, 0.02), ncol=1, framealpha = 1.0, numpoints = 1)

ax.text(215, -22.5, '$T^{\mathrm{max}} < 0^{\circ}$C  : '+str(ice_days_ref)+' yr$^{-1}$', color = 'dodgerblue', fontsize = 11, horizontalalignment='left', verticalalignment='center', zorder = 10)
ax.text(215, -27.5, '$T^{\mathrm{max}} \geq 20^{\circ}$C: '+str(warm_days_ref)+' yr$^{-1}$', color = 'yellowgreen', fontsize = 11, horizontalalignment='left', verticalalignment='center', zorder = 10)
ax.text(215, -32.5, '$T^{\mathrm{max}} \geq 25^{\circ}$C: '+str(summer_days_ref)+' yr$^{-1}$', color = 'sandybrown', fontsize = 11, horizontalalignment='left', verticalalignment='center', zorder = 10)
ax.text(215, -37.5, '$T^{\mathrm{max}} \geq 30^{\circ}$C: '+str(tropical_days_ref)+' yr$^{-1}$', color = 'indianred', fontsize = 11, horizontalalignment='left', verticalalignment='center', zorder = 10)

show()

