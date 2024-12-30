#Program plots the Arctic sea-ice fraction

from pylab import *
import numpy
import time
import glob, os
import math
import netCDF4 as netcdf
from scipy import stats
from cartopy import crs as ccrs, feature as cfeature
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as mticker

#Making pathway to folder with all data
directory_0600_PI	= '../../../Data/CESM_0600_PI/'

def ConverterField(index_break, field):
    """Shifts field, where it starts at 0E and ends at 360E"""

    new_field	= ma.masked_all(shape(field))
    length_section	= len(field[0]) - index_break

    #Shift the first part
    new_field[:, :length_section] = field[:, index_break:]

    #Shift the last part
    new_field[:, length_section:] = field[:, :index_break] 

    return new_field

def LowCESMPlot(lon, lat, field):
    """Returns 4 array's to plot on a global projection"""

    #Left of pole
    lon[lon > 180]	= lon[lon > 180] - 360.0

    lon_1		= lon[:, :160]
    lat_1		= lat[:, :160]
    field_1		= field[:, :160]

    #Right of pole
    lon_2		= lon[:, 159:]
    lat_2		= lat[:, 159:]
    field_2		= field[:, 159:]

    lat_3		= ma.masked_where(lon_2 > 0.0, lat_2)
    field_3		= ma.masked_where(lon_2 > 0.0, field_2)
    lon_3		= ma.masked_where(lon_2 > 0.0, lon_2)

    lon_2[lon_2 < -160] 	= lon_2[lon_2 < -160] + 360
    lat_2			= ma.masked_where(lon_2 < 0.0, lat_2)
    field_2			= ma.masked_where(lon_2 < 0.0, field_2)
    lon_2			= ma.masked_where(lon_2 < 0.0, lon_2)

    #To match at 40W
    index_1		= (fabs(lon[40] - 0.0)).argmin()

    lon_4		= ConverterField(index_1, lon)
    lat_4		= ConverterField(index_1, lat)
    field_4		= ConverterField(index_1, field)

    lon_4		= lon_4[:, 280:300]
    lat_4		= lat_4[:, 280:300]
    field_4		= field_4[:, 280:300]

    return lon_1, lat_1, field_1, lon_2, lat_2, field_2, lon_3, lat_3, field_3, lon_4, lat_4, field_4

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

year_start	= 1000
year_end	= 1099


#-----------------------------------------------------------------------------------------

fh 		= netcdf.Dataset(directory_0600_PI+'Ice/Arctic_ice_fraction_month_'+str(12)+'_year_'+str(year_start)+'-'+str(year_end)+'.nc', 'r')

lon		= fh.variables['lon'][:]
lat		= fh.variables['lat'][:] 
fraction_dec	= np.mean(fh.variables['Fraction'][:], axis = 0)

fh.close()

fh 		= netcdf.Dataset(directory_0600_PI+'Ice/Arctic_ice_fraction_month_'+str(1)+'_year_'+str(year_start)+'-'+str(year_end)+'.nc', 'r')

fraction_jan	= np.mean(fh.variables['Fraction'][:], axis = 0)

fh.close()

fh 		= netcdf.Dataset(directory_0600_PI+'Ice/Arctic_ice_fraction_month_'+str(2)+'_year_'+str(year_start)+'-'+str(year_end)+'.nc', 'r')

fraction_feb	= np.mean(fh.variables['Fraction'][:], axis = 0)

fh.close()

fh 		= netcdf.Dataset(directory_0600_PI+'Ice/Arctic_ice_fraction_month_'+str(3)+'_year_'+str(year_start)+'-'+str(year_end)+'.nc', 'r')

fraction_mar	= np.mean(fh.variables['Fraction'][:], axis = 0)

fh.close()

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

lon_1, lat_1, fraction_1_plot, lon_2, lat_2, fraction_2_plot, lon_3, lat_3, fraction_3_plot, lon_4, lat_4, fraction_4_plot	= LowCESMPlot(lon, lat, fraction_feb)

fig, ax = plt.subplots(subplot_kw={'projection': ccrs.Stereographic(central_longitude=-20, central_latitude=50)})

ax.fill_between([-80, 80], y1 = np.zeros(2) + 0, y2 = np.zeros(2) + 80, color = 'gray', alpha = 0.20, transform=ccrs.PlateCarree())

CS      = ax.contourf(lon_1, lat_1, fraction_1_plot, levels = np.arange(15, 100.01, 5), cmap = 'Blues_r', transform=ccrs.PlateCarree())
CS      = ax.contourf(lon_2, lat_2, fraction_2_plot, levels = np.arange(15, 100.01, 5), cmap = 'Blues_r', transform=ccrs.PlateCarree())
CS      = ax.contourf(lon_3, lat_3, fraction_3_plot, levels = np.arange(15, 100.01, 5), cmap = 'Blues_r', transform=ccrs.PlateCarree())
CS      = ax.contourf(lon_4, lat_4, fraction_4_plot, levels = np.arange(15, 100.01, 5), cmap = 'Blues_r', transform=ccrs.PlateCarree())

lon_1, lat_1, fraction_1_plot, lon_2, lat_2, fraction_2_plot, lon_3, lat_3, fraction_3_plot, lon_4, lat_4, fraction_4_plot	= LowCESMPlot(lon, lat, fraction_dec)

CS_1	= ax.contour(lon_1, lat_1, fraction_1_plot, levels = [15], colors = 'k', linewidths = 2, zorder = 20, transform=ccrs.PlateCarree())
CS_1	= ax.contour(lon_2, lat_2, fraction_2_plot, levels = [15], colors = 'k', linewidths = 2, zorder = 20,  transform=ccrs.PlateCarree())
CS_1	= ax.contour(lon_3, lat_3, fraction_3_plot, levels = [15], colors = 'k', linewidths = 2, zorder = 20,  transform=ccrs.PlateCarree())
CS_1	= ax.contour(lon_4, lat_4, fraction_4_plot, levels = [15], colors = 'k', linewidths = 2, zorder = 20,  transform=ccrs.PlateCarree())

lon_1, lat_1, fraction_1_plot, lon_2, lat_2, fraction_2_plot, lon_3, lat_3, fraction_3_plot, lon_4, lat_4, fraction_4_plot	= LowCESMPlot(lon, lat, fraction_jan)

CS_1	= ax.contour(lon_1, lat_1, fraction_1_plot, levels = [15], colors = 'firebrick', linewidths = 2, zorder = 20,  transform=ccrs.PlateCarree())
CS_1	= ax.contour(lon_2, lat_2, fraction_2_plot, levels = [15], colors = 'firebrick', linewidths = 2, zorder = 20,  transform=ccrs.PlateCarree())
CS_1	= ax.contour(lon_3, lat_3, fraction_3_plot, levels = [15], colors = 'firebrick', linewidths = 2, zorder = 20,  transform=ccrs.PlateCarree())
CS_1	= ax.contour(lon_4, lat_4, fraction_4_plot, levels = [15], colors = 'firebrick', linewidths = 2, zorder = 20,  transform=ccrs.PlateCarree())

divider = make_axes_locatable(ax)
ax_cb   = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
fig.add_axes(ax_cb)

cbar    = colorbar(CS, ticks = [15, 40, 60, 80, 100], cax=ax_cb)
cbar.set_label('Sea-ice fraction ($\%$)', fontsize = 12)

ax.add_feature(cfeature.LAND, zorder=1)
ax.coastlines()
ax.set_extent([-61, 11, 45, 78], ccrs.PlateCarree())

gl1	     = ax.gridlines(draw_labels=True, dms = True, x_inline=False, y_inline=False, linewidth = 0.0)
gl1.top_labels = False
gl1.right_labels = False
gl1.xlocator = mticker.FixedLocator([-40, -20, 0])
gl1.ylocator = mticker.FixedLocator([40, 50, 60])
gl1.xlabel_style = {'rotation':0}
gl2 	= ax.gridlines(draw_labels=False, dms = True, x_inline=False, y_inline=False)

graph_1	= ax.plot([50, 50], [90, 90], '-', color = 'k', linewidth = 2.0, label = 'December', transform=ccrs.PlateCarree(), zorder =10)
graph_2	= ax.plot([50, 50], [90, 90], '-', color = 'firebrick', linewidth = 2.0, label = 'January', transform=ccrs.PlateCarree(), zorder =10)

graphs	      = graph_1 + graph_2

legend_labels = [l.get_label() for l in graphs]
ax.legend(graphs, legend_labels, loc='lower right', ncol=1, framealpha = 1.0, numpoints = 1).set_zorder(12)

ax.set_title('a) Arctic sea-ice extent, February, AMOC on ($F_H$ = 0.18 Sv)')

show()

