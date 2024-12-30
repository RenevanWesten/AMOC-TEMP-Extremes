#Program plots the European January temperatures and 1:10-year cold extremes

from pylab import *
import numpy
import time
import glob, os
import math
import netCDF4 as netcdf
import matplotlib.colors as colors
from cartopy import crs as ccrs, feature as cfeature
from mpl_toolkits.axes_grid1 import make_axes_locatable

#Making pathway to folder with all data
directory_0600_PI	= '../../../Data/CESM_0600_PI/'
directory_3800_PI	= '../../../Data/CESM_3800_PI/'

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
    
    return lon_new, field_new
	
def ReturnValue(prob, shape, loc = 0.0, scale = 1.0):
    """Return the return value at a given probability"""
	
    return loc - (scale / shape)* (1.0 - (-np.log(1 - prob))**(-shape))

def ReturnTime(value, shape, loc = 0.0, scale = 1.0):
    """Returns the return time of a given event"""

    prob	= 1.0 - np.exp(-(1.0 + shape * ( (value - loc) / scale))**(-1.0 / shape))

    return 1.0 / prob
	
#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

month_start	    = 1
month_end	    = 1

return_period	= 10

#-----------------------------------------------------------------------------------------

fh = netcdf.Dataset(directory_3800_PI+'Atmosphere/TEMP_2m_Europe_month_'+str(month_start)+'-'+str(month_end)+'.nc', 'r')

time_all	= fh.variables['time'][:]
lon		    = fh.variables['lon'][:] 			
lat		    = fh.variables['lat'][:] 			
temp_all	= fh.variables['TEMP_2m'][:] 			 			

fh.close()

#-----------------------------------------------------------------------------------------

fh = netcdf.Dataset(directory_0600_PI+'Atmosphere/TEMP_2m_Europe_month_'+str(month_start)+'-'+str(month_end)+'.nc', 'r')

temp_ref_all	= fh.variables['TEMP_2m'][:] 			 			

fh.close()

#-----------------------------------------------------------------------------------------

lon_global, lat_global, return_level_min	= ReadinDataGEV(directory_3800_PI+'/Atmosphere/TEMP_2m_extremes_GEV_fit_minima.nc', return_period)

#-----------------------------------------------------------------------------------------
#Rescale the temperature plot
scale	= 5.0
cut_off	= 5

temp_plot			            = np.mean(temp_all, axis = 0) - np.mean(temp_ref_all, axis = 0)
temp_plot[temp_plot < -cut_off]	= (temp_plot[temp_plot < -cut_off] - -cut_off) / scale - cut_off
temp_plot[temp_plot > cut_off]	= (temp_plot[temp_plot > cut_off] - cut_off) / scale + cut_off

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})

CS      = ax.contourf(lon, lat, temp_plot, levels = np.arange(-12, 12.01, 0.5), extend = 'both', cmap = 'RdBu_r', transform=ccrs.PlateCarree())

divider = make_axes_locatable(ax)
ax_cb   = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
fig.add_axes(ax_cb)

cbar    = colorbar(CS, ticks = [-12, -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 12], cax=ax_cb)
cbar.ax.set_yticklabels([-40, -30, -20, -10, -4, -2, 0, 2, 4, 10, 20, 30, 40])
cbar.set_label('Temperature difference ($^{\circ}$C)')

CS_1	= ax.contour(lon_global, lat_global, return_level_min, levels = [-10], colors = 'firebrick', linewidths = 2, transform=ccrs.PlateCarree(), zorder = 10)
CS_1	= ax.contour(lon_global, lat_global, return_level_min, levels = [-20], colors = 'royalblue', linewidths = 2, transform=ccrs.PlateCarree(), zorder = 10)
CS_1	= ax.contour(lon_global, lat_global, -return_level_min, levels = [30], colors = 'k', linewidths = 2, transform=ccrs.PlateCarree(), zorder = 10)
CS_1	= ax.contour(lon_global, lat_global, -return_level_min, levels = [40], colors = 'b', linewidths = 2, transform=ccrs.PlateCarree(), zorder = 10)
CS_1	= ax.contour(lon_global, lat_global, -return_level_min, levels = [50], colors = 'r', linewidths = 2, transform=ccrs.PlateCarree(), zorder = 10)

gl = ax.gridlines(draw_labels=True)
gl.top_labels = False
gl.right_labels = False
ax.set_extent([-45, 45, 25, 75], ccrs.PlateCarree())
ax.coastlines('50m')
ax.add_feature(cfeature.LAND, zorder=0)

graph_1		= ax.plot([-100, -100], [-100, -100], '-', color = 'r', linewidth = 2, label = '$-50^{\circ}$C')
graph_2		= ax.plot([-100, -100], [-100, -100], '-', color = 'b', linewidth = 2, label = '$-40^{\circ}$C')
graph_3		= ax.plot([-100, -100], [-100, -100], '-', color = 'k', linewidth = 2, label = '$-30^{\circ}$C')
graph_5		= ax.plot([-100, -100], [-100, -100], '--', color = 'firebrick', linewidth = 2, label = '$-10^{\circ}$C')
graph_4		= ax.plot([-100, -100], [-100, -100], '--', color = 'royalblue', linewidth = 2, label = '$-20^{\circ}$C')

graphs	      	= graph_1 + graph_2 + graph_3 + graph_4 + graph_5
legend_labels 	= [l.get_label() for l in graphs]
legend_1	= ax.legend(graphs, legend_labels, loc='lower left', ncol=1, framealpha = 1.0, numpoints = 1).set_zorder(11)

ax.set_title('d) January temperature and 1:10-year cold extreme, AMOC off ($F_H = 0.18$ Sv)', fontsize = 11)

show()


