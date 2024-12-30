#Program plots the storm tracks

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
directory_1500_PI	= '../../../Data/CESM_1500_PI/'

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

month_start	= 1
month_end	= 1

#-----------------------------------------------------------------------------------------

fh = netcdf.Dataset(directory_1500_PI+'Atmosphere/TEMP_PSL_storm_tracks_month_'+str(month_start)+'-'+str(month_end)+'.nc', 'r')

lon         = fh.variables['lon'][:] 			
lat         = fh.variables['lat'][:] 
temp_var    = fh.variables['TEMP_var'][:]		
pres_var    = fh.variables['PSL_var'][:]

fh.close()

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})

CS      = ax.contourf(lon, lat, pres_var, levels = np.arange(0, 60.01, 5), extend = 'max', cmap = 'Spectral_r', transform=ccrs.PlateCarree())

ax.contour(lon, lat, temp_var, levels = [5], colors = 'b', linewidths = 2, zorder = 10)
ax.contour(lon, lat, temp_var, levels = [15], colors = 'k', linewidths = 2, zorder = 10)
ax.contour(lon, lat, temp_var, levels = [25], colors = 'r', linewidths = 2, zorder = 10)

divider = make_axes_locatable(ax)
ax_cb   = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
fig.add_axes(ax_cb)

cbar    = colorbar(CS, ticks = np.arange(0, 60.01, 10), cax=ax_cb)
cbar.set_label('Pressure variance (hPa$^{2}$)')

gl = ax.gridlines(draw_labels=True)
gl.top_labels = False
gl.right_labels = False
ax.set_extent([-90, 30, 19.999, 80.001], ccrs.PlateCarree())
ax.coastlines('110m')
ax.add_feature(cfeature.LAND, zorder=0)

graph_1	= ax.plot([50, 50], [90, 90], '-', color = 'r', linewidth = 2.0, label = '25 K$^2$', transform=ccrs.PlateCarree(), zorder = 10)
graph_2	= ax.plot([50, 50], [90, 90], '-', color = 'k', linewidth = 2.0, label = '15 K$^2$', transform=ccrs.PlateCarree(), zorder = 10)
graph_3	= ax.plot([50, 50], [90, 90], '-', color = 'b', linewidth = 2.0, label = '5 K$^2$', transform=ccrs.PlateCarree(), zorder = 10)

graphs	      = graph_1 + graph_2 + graph_3

legend_labels = [l.get_label() for l in graphs]
ax.legend(graphs, legend_labels, loc='lower left', ncol=1, framealpha = 1.0, numpoints = 1, prop={'size': 8}).set_zorder(12)

fig.subplots_adjust(top=0.969,bottom=0.031,left=0.086,right=0.859,hspace=0.2,wspace=0.2)

ax.set_title('c) Sea-level pressure and temperature variance, AMOC on ($F_H = 0.45$ Sv)', fontsize = 9)

show()
