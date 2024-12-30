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
import matplotlib.tri as tri
from scipy.interpolate import interp1d

#Making pathway to folder with all data
directory_1500_PI	= '../../../Data/CESM_1500_PI/'
directory_1500_RCP	= '../../../Data/CESM_1500_RCP45/'

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

month_start	= 1
month_end	= 1

#-----------------------------------------------------------------------------------------

fh = netcdf.Dataset(directory_1500_RCP+'Atmosphere/TEMP_PSL_storm_tracks_month_'+str(month_start)+'-'+str(month_end)+'.nc', 'r')

lon         = fh.variables['lon'][:] 			
lat         = fh.variables['lat'][:] 
temp_var    = fh.variables['TEMP_var'][:]		
pres_var    = fh.variables['PSL_var'][:]

fh.close()

#-----------------------------------------------------------------------------------------

fh = netcdf.Dataset(directory_1500_PI+'Atmosphere/TEMP_PSL_storm_tracks_month_'+str(month_start)+'-'+str(month_end)+'.nc', 'r')

temp_var_ref = fh.variables['TEMP_var'][:]
pres_var_ref = fh.variables['PSL_var'][:]

fh.close()

#-----------------------------------------------------------------------------------------

lon_loc, lat_loc	= 5.18, 52.10	#De Bilt

x, y                = np.meshgrid(lon, lat)
x, y                = x.ravel(), y.ravel()
triang              = tri.Triangulation(x, y)

temp_int            = tri.LinearTriInterpolator(triang, temp_var.ravel())
pres_int            = tri.LinearTriInterpolator(triang, pres_var.ravel())
temp_ref_int        = tri.LinearTriInterpolator(triang, temp_var_ref.ravel())
pres_ref_int        = tri.LinearTriInterpolator(triang, pres_var_ref.ravel())

print('2 to 6-day variability for De Bilt:')
print('Temperature                 : ', temp_int(lon_loc, lat_loc), 'K^2')
print('Sea-level pressure          : ', pres_int(lon_loc, lat_loc), 'hPa^2')
print()
print('Temperature (AMOC on)       : ', temp_ref_int(lon_loc, lat_loc), 'K^2')
print('Sea-level pressure (AMOC on): ', pres_ref_int(lon_loc, lat_loc), 'hPa^2')

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})

CS      = ax.contourf(lon, lat, pres_var - pres_var_ref, levels = np.arange(-50, 50.1, 5), extend = 'both', cmap = 'RdBu_r', transform=ccrs.PlateCarree())

ax.contour(lon, lat, (temp_var - temp_var_ref), levels = [-20], colors = 'royalblue', linewidths = 2, zorder = 10)
ax.contour(lon, lat, -(temp_var - temp_var_ref), levels = [10], colors = 'royalblue', linewidths = 2, zorder = 10)
ax.contour(lon, lat, temp_var - temp_var_ref, levels = [0], colors = 'k', linewidths = 2, zorder = 10)
ax.contour(lon, lat, temp_var - temp_var_ref, levels = [10], colors = 'firebrick', linewidths = 2, zorder = 10)
ax.contour(lon, lat, -(temp_var - temp_var_ref), levels = [-20], colors = 'firebrick', linewidths = 2, zorder = 10)

divider = make_axes_locatable(ax)
ax_cb   = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
fig.add_axes(ax_cb)

cbar    = colorbar(CS, ticks = np.arange(-50, 50.01, 10), cax=ax_cb)
cbar.set_label('Pressure variance difference (hPa$^{2}$)')

gl = ax.gridlines(draw_labels=True)
gl.top_labels = False
gl.right_labels = False
ax.set_extent([-90, 30, 19.999, 80.001], ccrs.PlateCarree())
ax.coastlines('110m')
ax.add_feature(cfeature.LAND, zorder=0)

graph_1	= ax.plot([50, 50], [90, 90], '--', color = 'firebrick', linewidth = 2.0, label = '$+$20 K$^2$', transform=ccrs.PlateCarree(), zorder = 10)
graph_2	= ax.plot([50, 50], [90, 90], '-', color = 'firebrick', linewidth = 2.0, label = '$+$10 K$^2$', transform=ccrs.PlateCarree(), zorder = 10)
graph_3	= ax.plot([50, 50], [90, 90], '-', color = 'k', linewidth = 2.0, label = '0 K$^2$', transform=ccrs.PlateCarree(), zorder = 10)
graph_4	= ax.plot([50, 50], [90, 90], '-', color = 'royalblue', linewidth = 2.0, label = '$-$10 K$^2$', transform=ccrs.PlateCarree(), zorder = 10)
graph_5	= ax.plot([50, 50], [90, 90], '--', color = 'royalblue', linewidth = 2.0, label = '$-$20 K$^2$', transform=ccrs.PlateCarree(), zorder = 10)

graphs	      = graph_1 + graph_2 + graph_3 + graph_4 + graph_5

legend_labels = [l.get_label() for l in graphs]
ax.legend(graphs, legend_labels, loc='lower left', ncol=1, framealpha = 1.0, numpoints = 1, prop={'size': 8}).set_zorder(12)

fig.subplots_adjust(top=0.969,bottom=0.031,left=0.086,right=0.859,hspace=0.2,wspace=0.2)

ax.set_title('k) Sea-level pressure and temperature variance, RCP4.5 ($F_H = 0.45$ Sv)', fontsize = 9)

show()

