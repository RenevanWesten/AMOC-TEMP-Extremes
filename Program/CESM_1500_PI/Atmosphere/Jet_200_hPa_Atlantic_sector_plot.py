#Program plots the wind speeds at 200 hPa

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

fh = netcdf.Dataset(directory_1500_PI+'Atmosphere/Jet_200_hPa_Atlantic_sector_month_'+str(month_start)+'-'+str(month_end)+'.nc', 'r')

time_all	= fh.variables['time'][:]
lon		    = fh.variables['lon'][:] 			
lat		    = fh.variables['lat'][:] 			
u_vel_all	= fh.variables['U'][:]
u_vel_2_all	= fh.variables['UU'][:] 			
v_vel_all	= fh.variables['V'][:] 				
v_vel_2_all	= fh.variables['VV'][:]	

fh.close()

#Get the wind speed mean and direction
vel_speed	= np.mean(np.sqrt(u_vel_2_all + v_vel_2_all), axis = 0)
u_vel_all	= np.mean(u_vel_all, axis = 0)
v_vel_all	= np.mean(v_vel_all, axis = 0)

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})

CS      = ax.contourf(lon, lat, vel_speed, levels = np.arange(0, 50.01, 2.5), extend = 'max', cmap = 'Spectral_r', transform=ccrs.PlateCarree())

divider = make_axes_locatable(ax)
ax_cb   = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
fig.add_axes(ax_cb)

cbar    = colorbar(CS, ticks = np.arange(0, 50.01, 10), cax=ax_cb)
cbar.set_label('Wind speed (m s$^{-1}$)')

scale_arrow	= 4
Q = ax.quiver(lon[::scale_arrow], lat[::scale_arrow], u_vel_all[::scale_arrow, ::scale_arrow], v_vel_all[::scale_arrow, ::scale_arrow], scale = 600, transform=ccrs.PlateCarree())

qk = ax.quiverkey(Q, 0.17, 0.10, 40, '40 m s$^{-1}$', labelpos = 'S', coordinates='figure')

ax.plot([-45, -45], [-10, 85], '--', linewidth = 2.0, color = 'royalblue')
ax.plot([15, 15], [-10, 85], '--', linewidth = 2.0, color = 'royalblue')

gl = ax.gridlines(draw_labels=True)
gl.top_labels = False
gl.right_labels = False
ax.set_extent([-90, 30, -0.001, 80.001], ccrs.PlateCarree())
ax.coastlines('110m')
ax.add_feature(cfeature.LAND, zorder=0)

ax.set_title('c) 200 hPa velocities, AMOC on ($F_H = 0.45$ Sv)')

show()


