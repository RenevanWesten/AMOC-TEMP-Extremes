#Program plots the 1:10-year temperature extremes

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
import matplotlib.colors as colors
from cartopy import crs as ccrs, feature as cfeature
from mpl_toolkits.axes_grid1 import make_axes_locatable

#Making pathway to folder with all data
directory_1500_PI	= '../../../Data/CESM_1500_PI/'
directory_1500_RCP	= '../../../Data/CESM_1500_RCP45/'

def ReadinDataGEV(filename, return_period):

    fh = netcdf.Dataset(filename, 'r')

    lon		    = fh.variables['lon'][:]
    lat		    = fh.variables['lat'][:]
    shape_all	= fh.variables['shape'][:] 
    loc_all		= fh.variables['loc'][:] 
    scale_all	= fh.variables['scale'][:] 

    fh.close()

    #-----------------------------------------------------------------------------------------

    return_level	= ma.masked_all((len(lat), len(lon)))

    for lat_i in range(len(lat)):
        for lon_i in range(len(lon)):
            #Determine the return level the data for the fourth fit
            return_level[lat_i, lon_i]     = ReturnValue(1.0/return_period, shape_all[lat_i, lon_i], loc_all[lat_i, lon_i], scale_all[lat_i, lon_i]) #Fit curve

    lon, return_level	= ConverterField(lon, lat, return_level)

    return lon, lat, return_level

def ConverterField(lon, lat, field):
    """Shifts field, to -180E to 180E"""
    lon_new		= ma.masked_all(len(lon))
    field_new	= ma.masked_all(np.shape(field))

    #Get the corresponding index
    index		= (fabs(lon - 180.0)).argmin()

    #Start filling at -180E
    lon_new[:len(lon[index:])] 	 	    = lon[index:] - 360.0
    field_new[:, :len(lon[index:])]	  	= field[:, index:]

    #Fill the remaining part to 360E
    lon_new[len(lon[index:]):] 		    = lon[:index]
    field_new[:, len(lon[index:]):]		= field[:, :index]
    
    lon_new, field_new	= PeriodicBoundaries2D(lon_new, lat, field_new)
    
    return lon_new, field_new
	
def ReturnValue(prob, shape, loc = 0.0, scale = 1.0):
    """Return the return value at a given probability"""
	
    return loc - (scale / shape)* (1.0 - (-np.log(1 - prob))**(-shape))

def ReturnTime(value, shape, loc = 0.0, scale = 1.0):
    """Returns the return time of a given event"""

    prob	= 1.0 - np.exp(-(1.0 + shape * ( (value - loc) / scale))**(-1.0 / shape))

    return 1.0 / prob

def PeriodicBoundaries2D(lon, lat, field, lon_grids = 1):
    """Add periodic zonal boundaries for 2D field"""

    #Empty field with additional zonal boundaries
    lon_2			= np.zeros(len(lon) + lon_grids * 2)
    field_2			= ma.masked_all((len(lat), len(lon_2)))
	
    #Get the left boundary, which is the right boundary of the original field
    lon_2[:lon_grids]	= lon[-lon_grids:] - 360.0
    field_2[:, :lon_grids]	= field[:, -lon_grids:]

    #Same for the right boundary
    lon_2[-lon_grids:]	= lon[:lon_grids] + 360.0
    field_2[:, -lon_grids:]	= field[:, :lon_grids]

    #And the complete field
    lon_2[lon_grids:-lon_grids]		= lon
    field_2[:, lon_grids:-lon_grids] 	= field

    return lon_2, field_2	


#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

return_period		= 10	#In years

#-----------------------------------------------------------------------------------------

lon, lat, return_level_min_1500	= ReadinDataGEV(directory_1500_PI+'/Atmosphere/TEMP_2m_extremes_GEV_fit_minima.nc', return_period)
lon, lat, return_level_max_1500	= ReadinDataGEV(directory_1500_PI+'/Atmosphere/TEMP_2m_extremes_GEV_fit_maxima.nc', return_period)

lon, lat, return_level_min_RCP	= ReadinDataGEV(directory_1500_RCP+'/Atmosphere/TEMP_2m_extremes_GEV_fit_minima.nc', return_period)
lon, lat, return_level_max_RCP	= ReadinDataGEV(directory_1500_RCP+'/Atmosphere/TEMP_2m_extremes_GEV_fit_maxima.nc', return_period)

#-----------------------------------------------------------------------------------------
#Rescale the temperature plot
scale	= 5.0
cut_off	= 5

return_level_min_plot					= return_level_min_RCP - return_level_min_1500
return_level_min_plot[return_level_min_plot < -cut_off]	= (return_level_min_plot[return_level_min_plot < -cut_off] - -cut_off) / scale - cut_off
return_level_min_plot[return_level_min_plot > cut_off]	= (return_level_min_plot[return_level_min_plot > cut_off] - cut_off) / scale + cut_off

return_level_max_plot					= return_level_max_RCP - return_level_max_1500
return_level_max_plot[return_level_max_plot < -cut_off]	= (return_level_max_plot[return_level_max_plot < -cut_off] - -cut_off) / scale - cut_off
return_level_max_plot[return_level_max_plot > cut_off]	= (return_level_max_plot[return_level_max_plot > cut_off] - cut_off) / scale + cut_off

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax = plt.subplots(subplot_kw={'projection': ccrs.Robinson()})

CS      = ax.contourf(lon, lat, return_level_min_plot, levels = np.arange(-12, 12.01, 0.5), extend = 'both', cmap = 'RdBu_r', transform=ccrs.PlateCarree())

divider = make_axes_locatable(ax)
ax_cb   = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
fig.add_axes(ax_cb)

cbar    = colorbar(CS, ticks = [-12, -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 12], cax=ax_cb)
cbar.ax.set_yticklabels([-40, -30, -20, -10, -4, -2, 0, 2, 4, 10, 20, 30, 40])
cbar.set_label('Temperature difference ($^{\circ}$C)')

ax.set_global()
ax.gridlines()
ax.coastlines()

for lat_i in range(0, len(lat), 3):
    for lon_i in range(0, len(lon), 3):

        if return_level_min_RCP[lat_i, lon_i] <= -40:
            #Non-significant difference
            ax.scatter(lon[lon_i], lat[lat_i], marker = 'o', edgecolor = 'k' , s = 6, facecolors='none', transform=ccrs.PlateCarree())

ax.set_title('k) 1:10-year cold extreme, RCP4.5 ($F_H = 0.45$ Sv)')

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax = plt.subplots(subplot_kw={'projection': ccrs.Robinson()})

CS      = ax.contourf(lon, lat, return_level_max_plot, levels = np.arange(-12, 12.01, 0.5), extend = 'both', cmap = 'RdBu_r', transform=ccrs.PlateCarree())

divider = make_axes_locatable(ax)
ax_cb   = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
fig.add_axes(ax_cb)

cbar    = colorbar(CS, ticks = [-12, -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 12], cax=ax_cb)
cbar.ax.set_yticklabels([-40, -30, -20, -10, -4, -2, 0, 2, 4, 10, 20, 30, 40])
cbar.set_label('Temperature difference ($^{\circ}$C)')

ax.set_global()
ax.gridlines()
ax.coastlines()

for lat_i in range(0, len(lat), 3):
    for lon_i in range(0, len(lon), 3):

        if return_level_max_RCP[lat_i, lon_i] >= 40:
            #Non-significant difference
            ax.scatter(lon[lon_i], lat[lat_i], marker = 'o', edgecolor = 'k' , s = 6, facecolors='none', transform=ccrs.PlateCarree())

ax.set_title('l) 1:10-year warm extreme, RCP4.5 ($F_H = 0.45$ Sv)')

show()

