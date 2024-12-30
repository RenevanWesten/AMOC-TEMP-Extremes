#Program plots the atmospheric meridional heat transports

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
directory_1500_RCP	= '../../../Data/CESM_1500_RCP85/'

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------	

month_start	= 1	
month_end	= 1	

#-----------------------------------------------------------------------------------------	
#-----------------------------------------------------------------------------------------	
#-----------------------------------------------------------------------------------------	
	
fh = netcdf.Dataset(directory_1500_RCP+'Atmosphere/Meridional_SHF_LHF_month_'+str(month_start)+'-'+str(month_end)+'.nc', 'r')

lat             = fh.variables['lat'][:]     		
SHF_global      = fh.variables['SHF'][:]     		
SHF_eddy_global = fh.variables['SHF_eddy'][:]   
LHF_global      = fh.variables['LHF'][:]     		
LHF_eddy_global	= fh.variables['LHF_eddy'][:]     	
SHF_atl         = fh.variables['SHF_ATL'][:]     
SHF_eddy_atl	= fh.variables['SHF_eddy_ATL'][:]  	
LHF_atl         = fh.variables['LHF_ATL'][:]     		
LHF_eddy_atl    = fh.variables['LHF_eddy_ATL'][:]  	

fh.close()

#-----------------------------------------------------------------------------------------

fh = netcdf.Dataset(directory_1500_PI+'Atmosphere/Meridional_SHF_LHF_month_'+str(month_start)+'-'+str(month_end)+'.nc', 'r')

SHF_global_ref		= fh.variables['SHF'][:]     		
SHF_eddy_global_ref	= fh.variables['SHF_eddy'][:]   
LHF_global_ref		= fh.variables['LHF'][:]     		
LHF_eddy_global_ref	= fh.variables['LHF_eddy'][:]     	
SHF_atl_ref		    = fh.variables['SHF_ATL'][:]     
SHF_eddy_atl_ref	= fh.variables['SHF_eddy_ATL'][:]  	
LHF_atl_ref		    = fh.variables['LHF_ATL'][:]     		
LHF_eddy_atl_ref	= fh.variables['LHF_eddy_ATL'][:]  	

fh.close()

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax		= subplots()

ax.plot(lat, SHF_global - SHF_global_ref, '-', color = 'r', linewidth = 2.0)
ax.plot(lat, SHF_atl - SHF_atl_ref, '--', color = 'firebrick', linewidth = 2.0)

ax.plot(lat, LHF_global - LHF_global_ref, '-', color = 'b', linewidth = 2.0)
ax.plot(lat, LHF_atl - LHF_atl_ref, '--', color = 'royalblue', linewidth = 2.0)

ax.set_ylabel('Meridional heat transport difference (PW)')
ax.set_xlim(0, 80)
ax.set_ylim(-17, 17)
ax.grid()

ax.set_xticks(np.arange(0, 80.1, 20))
ax.set_xticklabels(['Eq', '20$^{\circ}$N', '40$^{\circ}$N', '60$^{\circ}$N', '80$^{\circ}$N'])

graph_1		= ax.plot([90, 90], [-1, -1], linestyle = '-', color = 'r', linewidth = 2, label = 'MHT$_{\mathrm{S}}$ (Global)')
graph_2		= ax.plot([90, 90], [-1, -1], linestyle = '--', color = 'firebrick', linewidth = 2, label = 'MHT$_{\mathrm{S}}$ (Atlantic)')
graph_3		= ax.plot([90, 90], [-1, -1], linestyle = '-', color = 'b', linewidth = 2, label = 'MHT$_{\mathrm{L}}$ (Global)')
graph_4		= ax.plot([90, 90], [-1, -1], linestyle = '--', color = 'royalblue', linewidth = 2, label = 'MHT$_{\mathrm{L}}$ (Atlantic)')

graphs		= graph_1 + graph_2 + graph_3 + graph_4
legend_labels = [l.get_label() for l in graphs]
legend_1      = ax.legend(graphs, legend_labels, loc = 'lower left', ncol=1, numpoints = 1, framealpha = 1.0).set_zorder(20)

ax.set_title('o) MHT (January), RCP8.5 ($F_H = 0.45$ Sv)')

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax		= subplots()

ax.plot(lat, SHF_eddy_global - SHF_eddy_global_ref, '-', color = 'r', linewidth = 2.0)
ax.plot(lat, SHF_eddy_atl - SHF_eddy_atl_ref, '--', color = 'firebrick', linewidth = 2.0)

ax.plot(lat, LHF_eddy_global - LHF_eddy_global_ref, '-', color = 'b', linewidth = 2.0)
ax.plot(lat, LHF_eddy_atl - LHF_eddy_atl_ref, '--', color = 'royalblue', linewidth = 2.0)

ax.set_ylabel('Meridional heat transport difference (PW)')
ax.set_xlim(0, 80)
ax.set_ylim(-1, 1)
ax.set_yticks([-1, -0.5, 0, 0.5, 1])
ax.grid()

ax.set_xticks(np.arange(0, 80.1, 20))
ax.set_xticklabels(['Eq', '20$^{\circ}$N', '40$^{\circ}$N', '60$^{\circ}$N', '80$^{\circ}$N'])

graph_1		= ax.plot([90, 90], [-1, -1], linestyle = '-', color = 'r', linewidth = 2, label = 'MHT$_{\mathrm{S,eddy}}$ (Global)')
graph_2		= ax.plot([90, 90], [-1, -1], linestyle = '--', color = 'firebrick', linewidth = 2, label = 'MHT$_{\mathrm{S,eddy}}$ (Atlantic)')
graph_3		= ax.plot([90, 90], [-1, -1], linestyle = '-', color = 'b', linewidth = 2, label = 'MHT$_{\mathrm{L,eddy}}$ (Global)')
graph_4		= ax.plot([90, 90], [-1, -1], linestyle = '--', color = 'royalblue', linewidth = 2, label = 'MHT$_{\mathrm{L,eddy}}$ (Atlantic)')

graphs		= graph_1 + graph_2 + graph_3 + graph_4
legend_labels = [l.get_label() for l in graphs]
legend_1      = ax.legend(graphs, legend_labels, loc = 'lower right', ncol=1, numpoints = 1, framealpha = 1.0).set_zorder(20)

ax.set_title('p) Eddy MHT (January), RCP8.5 ($F_H = 0.45$ Sv)')

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax		= subplots(figsize = (3.5, 4.5))

ax.plot(SHF_eddy_atl - SHF_eddy_atl_ref, lat, '--', color = 'firebrick', linewidth = 2.0)
ax.plot(LHF_eddy_atl - LHF_eddy_atl_ref, lat, '--', color = 'royalblue', linewidth = 2.0)

ax.set_xlabel('MHT difference (PW)', fontsize = 14)
ax.set_ylim(20, 80)
ax.set_xlim(-0.4, 0.4)
ax.grid()
plt.yticks(fontsize=14)
plt.xticks(fontsize=14)

ax.set_yticks(np.arange(20, 80.1, 10))
ax.set_yticklabels(['20$^{\circ}$N', '30$^{\circ}$N','40$^{\circ}$N','50$^{\circ}$N','60$^{\circ}$N','70$^{\circ}$N','80$^{\circ}$N'])

graph_1		= ax.plot([90, 90], [-1, -1], linestyle = '--', color = 'firebrick', linewidth = 2, label = 'MHT$_{\mathrm{S,eddy}}$')
graph_2		= ax.plot([90, 90], [-1, -1], linestyle = '--', color = 'royalblue', linewidth = 2, label = 'MHT$_{\mathrm{L,eddy}}$')

graphs		= graph_1 + graph_2
legend_labels = [l.get_label() for l in graphs]
legend_1      = ax.legend(graphs, legend_labels, loc = 'lower left', ncol=1, numpoints = 1, framealpha = 1.0).set_zorder(20)

fig.subplots_adjust(top=0.915,bottom=0.132,left=0.171,right=0.926,hspace=0.2,wspace=0.2)

fig.suptitle('p) Eddy MHT, RCP8.5 ($F_H = 0.45$ Sv)', fontsize = 13)

#-----------------------------------------------------------------------------------------	
#-----------------------------------------------------------------------------------------	
#-----------------------------------------------------------------------------------------	
	
fh = netcdf.Dataset(directory_1500_RCP+'Atmosphere/SHF_LHF_month_'+str(month_start)+'-'+str(month_end)+'.nc', 'r')

lon         = fh.variables['lon'][:]     		
lat         = fh.variables['lat'][:]     		
SHF         = fh.variables['SHF'][:]     		
SHF_eddy    = fh.variables['SHF_eddy'][:]   
LHF         = fh.variables['LHF'][:]     		
LHF_eddy    = fh.variables['LHF_eddy'][:]     	

fh.close()

#-----------------------------------------------------------------------------------------	
	
fh = netcdf.Dataset(directory_1500_PI+'Atmosphere/SHF_LHF_month_'+str(month_start)+'-'+str(month_end)+'.nc', 'r')
   		
SHF_ref		= fh.variables['SHF'][:]     		
SHF_eddy_ref= fh.variables['SHF_eddy'][:]   
LHF_ref		= fh.variables['LHF'][:]     		
LHF_eddy_ref= fh.variables['LHF_eddy'][:]     	

fh.close()

#-----------------------------------------------------------------------------------------

SHF_plot        = np.mean(SHF, axis = 0) - np.mean(SHF_ref, axis = 0)
SHF_eddy_plot   = np.mean(SHF_eddy, axis = 0) - np.mean(SHF_eddy_ref, axis = 0)
LHF_plot        = np.mean(LHF, axis = 0) - np.mean(LHF_ref, axis = 0)
LHF_eddy_plot   = np.mean(LHF_eddy, axis = 0) - np.mean(LHF_eddy_ref, axis = 0)

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})

CS      = ax.contourf(lon, lat, SHF_plot, levels = np.arange(-2, 2.01, 0.2), extend = 'both', cmap = 'RdBu_r', transform=ccrs.PlateCarree())

divider = make_axes_locatable(ax)
ax_cb   = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
fig.add_axes(ax_cb)

cbar    = colorbar(CS, ticks = np.arange(-2, 2.01, 1), cax=ax_cb)
cbar.set_label('Meridional heat transport difference (PW)')

ax.plot([-45, -45], [-10, 85], '--', linewidth = 2.0, color = 'royalblue')
ax.plot([15, 15], [-10, 85], '--', linewidth = 2.0, color = 'royalblue')

gl = ax.gridlines(draw_labels=True)
gl.top_labels = False
gl.right_labels = False
ax.set_extent([-90, 30, -0.001, 80.001], ccrs.PlateCarree())
ax.coastlines('110m')
ax.add_feature(cfeature.LAND, zorder=0)

fig.subplots_adjust(top=0.969,bottom=0.031,left=0.086,right=0.864,hspace=0.2,wspace=0.2)

ax.set_title('o) MHT$_{\mathrm{S}}$ (January), RCP8.5 ($F_H = 0.45$ Sv)')

#-----------------------------------------------------------------------------------------

fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})

CS      = ax.contourf(lon, lat, SHF_eddy_plot, levels = np.arange(-0.02, 0.021, 0.002), extend = 'both', cmap = 'RdBu_r', transform=ccrs.PlateCarree())

divider = make_axes_locatable(ax)
ax_cb   = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
fig.add_axes(ax_cb)

cbar    = colorbar(CS, ticks = np.arange(-0.02, 0.021, 0.01), cax=ax_cb)
cbar.set_label('Meridional heat transport difference (PW)')

ax.plot([-45, -45], [-10, 85], '--', linewidth = 2.0, color = 'royalblue')
ax.plot([15, 15], [-10, 85], '--', linewidth = 2.0, color = 'royalblue')

gl = ax.gridlines(draw_labels=True)
gl.top_labels = False
gl.right_labels = False
ax.set_extent([-90, 30, -0.001, 80.001], ccrs.PlateCarree())
ax.coastlines('110m')
ax.add_feature(cfeature.LAND, zorder=0)

fig.subplots_adjust(top=0.969,bottom=0.031,left=0.086,right=0.864,hspace=0.2,wspace=0.2)

ax.set_title('p) MHT$_{\mathrm{S,eddy}}$ (January), RCP8.5 ($F_H = 0.45$ Sv)')

#-----------------------------------------------------------------------------------------

fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})

CS      = ax.contourf(lon, lat, LHF_plot, levels = np.arange(-2, 2.01, 0.2), extend = 'both', cmap = 'RdBu_r', transform=ccrs.PlateCarree())

divider = make_axes_locatable(ax)
ax_cb   = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
fig.add_axes(ax_cb)

cbar    = colorbar(CS, ticks = np.arange(-2, 2.01, 1), cax=ax_cb)
cbar.set_label('Meridional heat transport difference (PW)')

ax.plot([-45, -45], [-10, 85], '--', linewidth = 2.0, color = 'royalblue')
ax.plot([15, 15], [-10, 85], '--', linewidth = 2.0, color = 'royalblue')

gl = ax.gridlines(draw_labels=True)
gl.top_labels = False
gl.right_labels = False
ax.set_extent([-90, 30, -0.001, 80.001], ccrs.PlateCarree())
ax.coastlines('110m')
ax.add_feature(cfeature.LAND, zorder=0)

fig.subplots_adjust(top=0.969,bottom=0.031,left=0.086,right=0.864,hspace=0.2,wspace=0.2)

ax.set_title('MHT by LHF (January), RCP8.5 ($F_H = 0.45$ Sv)')

#-----------------------------------------------------------------------------------------

fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})

CS      = ax.contourf(lon, lat, LHF_eddy_plot, levels = np.arange(-0.02, 0.021, 0.002), extend = 'both', cmap = 'RdBu_r', transform=ccrs.PlateCarree())

divider = make_axes_locatable(ax)
ax_cb   = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
fig.add_axes(ax_cb)

cbar    = colorbar(CS, ticks = np.arange(-0.02, 0.021, 0.01), cax=ax_cb)
cbar.set_label('Meridional heat transport difference (PW)')

ax.plot([-45, -45], [-10, 85], '--', linewidth = 2.0, color = 'royalblue')
ax.plot([15, 15], [-10, 85], '--', linewidth = 2.0, color = 'royalblue')

gl = ax.gridlines(draw_labels=True)
gl.top_labels = False
gl.right_labels = False
ax.set_extent([-90, 30, -0.001, 80.001], ccrs.PlateCarree())
ax.coastlines('110m')
ax.add_feature(cfeature.LAND, zorder=0)

fig.subplots_adjust(top=0.969,bottom=0.031,left=0.086,right=0.864,hspace=0.2,wspace=0.2)

ax.set_title('Eddy MHT by LHF (January), RCP8.5 ($F_H = 0.45$ Sv)')

show()
