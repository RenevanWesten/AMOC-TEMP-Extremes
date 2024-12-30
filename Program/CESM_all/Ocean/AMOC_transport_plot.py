#Program plots the AMOC strength at 1,000 m depth and 26N

from pylab import *
import numpy
import time
import glob, os
import math
import netCDF4 as netcdf
import matplotlib.colors as colors

#Making pathway to folder with all data
directory	  = '../../../Data/'

def ReadinData(filename):

    fh = netcdf.Dataset(filename, 'r')
    
    time		= fh.variables['time'][:]		
    transport	= fh.variables['Transport'][:]	#AMOC strength (Sv)
    
    fh.close()
    
    return time, transport
	
#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------	

depth_min 	= 0
depth_max	= 1000

#-----------------------------------------------------------------------------------------	

#The steady states under constant F_H and PI conditions
time_0600, transport_0600	= ReadinData(directory+'CESM_0600_PI/Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m.nc')	
time_1500, transport_1500	= ReadinData(directory+'CESM_1500_PI/Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m.nc')	
time_2900, transport_2900	= ReadinData(directory+'CESM_2900_PI/Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m.nc')	
time_3800, transport_3800	= ReadinData(directory+'CESM_3800_PI/Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m.nc')	

#The climate change simulations
time_hist_rcp45_branch_0600, transport_hist_rcp45_branch_0600		= ReadinData(directory+'CESM_0600_RCP45/Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m.nc')
time_hist_rcp85_branch_0600, transport_hist_rcp85_branch_0600		= ReadinData(directory+'CESM_0600_RCP85/Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m.nc')
time_hist_rcp45_branch_1500, transport_hist_rcp45_branch_1500		= ReadinData(directory+'CESM_1500_RCP45/Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m.nc')
time_hist_rcp85_branch_1500, transport_hist_rcp85_branch_1500		= ReadinData(directory+'CESM_1500_RCP85/Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m.nc')

#Model year 0 is from the quasi-equilibrium, model year 1 is the first year with constant F_H
time_branch	                = np.arange(len(transport_0600))

#The climate change simulations start at model year 500 from the steady states
time_branch_climate_change	= np.arange(500, 500+len(time_hist_rcp45_branch_0600))

#The quasi-equilbrium simulation
time_QE, transport_QE		= ReadinData(directory+'CESM_QE/Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m.nc')

#-----------------------------------------------------------------------------------------	
#-----------------------------------------------------------------------------------------	
#-----------------------------------------------------------------------------------------	

fig, ax	= subplots()

ax.fill_between([401, 500], -5, 25, alpha=0.25, edgecolor='orange', facecolor='orange')
ax.fill_between([1050, 1149], -5, 25, alpha=0.25, edgecolor='orange', facecolor='orange')

ax.plot(time_branch, transport_0600, '-', color = 'r', linewidth = 0.5)
ax.plot(time_branch, transport_3800, '-', color = 'b', linewidth = 0.5)

ax.text(450, 18, 'c', verticalalignment='center', horizontalalignment='center', color = 'r', fontsize=16, fontweight='bold')
ax.text(450, 7, 'd', verticalalignment='center', horizontalalignment='center', color = 'b', fontsize=16, fontweight='bold')
ax.text(1100, 18, 'e', verticalalignment='center', horizontalalignment='center', color = 'royalblue', fontsize=16, fontweight='bold')
ax.text(1100, 1, 'f', verticalalignment='center', horizontalalignment='center', color = 'firebrick', fontsize=16, fontweight='bold')

ax.plot(time_branch_climate_change[155:], transport_hist_rcp45_branch_0600[155:], '-', color = 'royalblue', linewidth = 0.5)
ax.plot(time_branch_climate_change[155:], transport_hist_rcp85_branch_0600[155:], '-', color = 'firebrick', linewidth = 0.5)
ax.plot(time_branch_climate_change[:156], transport_hist_rcp45_branch_0600[:156], '-', color = 'k', linewidth = 0.5)

ax.set_xlim(0, 1200)
ax.set_ylim(-2, 22)
ax.set_xlabel('Model year after branching point')
ax.set_ylabel('Volume transport (Sv)')
ax.grid()

graph_1		= ax.plot([-100, -100], [-100, -100], '-', color = 'r', linewidth = 2, label = 'AMOC on')
graph_2		= ax.plot([-100, -100], [-100, -100], '-', color = 'b', linewidth = 2, label = 'AMOC off')
graph_3		= ax.plot([-100, -100], [-100, -100], '-', color = 'k', linewidth = 2, label = 'Historical')
graph_4		= ax.plot([-100, -100], [-100, -100], '-', color = 'royalblue', linewidth = 2, label = 'RCP4.5')
graph_5		= ax.plot([-100, -100], [-100, -100], '-', color = 'firebrick', linewidth = 2, label = 'RCP8.5')

graphs	      	= graph_1 + graph_2 + graph_3 + graph_4 + graph_5
legend_labels 	= [l.get_label() for l in graphs]
legend_1	= ax.legend(graphs, legend_labels, loc=(0.76, 0.31), ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_title('a) AMOC strength at 26$^{\circ}$N ($F_H = 0.18$ Sv)')

#-----------------------------------------------------------------------------------------	
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.fill_between([401, 500], -5, 25, alpha=0.25, edgecolor='orange', facecolor='orange')
ax.fill_between([1050, 1149], -5, 25, alpha=0.25, edgecolor='orange', facecolor='orange')

ax.plot(time_branch, transport_1500, '-', color = 'r', linewidth = 0.5)
ax.plot(time_branch, transport_2900, '-', color = 'b', linewidth = 0.5)

ax.text(450, 14.5, 'g', verticalalignment='center', horizontalalignment='center', color = 'r', fontsize=16, fontweight='bold')
ax.text(450, 2.5, 'h', verticalalignment='center', horizontalalignment='center', color = 'b', fontsize=16, fontweight='bold')
ax.text(1100, 3, 'i', verticalalignment='center', horizontalalignment='center', color = 'royalblue', fontsize=16, fontweight='bold')
ax.text(1100, -1, 'j', verticalalignment='center', horizontalalignment='center', color = 'firebrick', fontsize=16, fontweight='bold')

ax.plot(time_branch_climate_change[155:], transport_hist_rcp45_branch_1500[155:], '-', color = 'royalblue', linewidth = 0.5)
ax.plot(time_branch_climate_change[155:], transport_hist_rcp85_branch_1500[155:], '-', color = 'firebrick', linewidth = 0.5)
ax.plot(time_branch_climate_change[:156], transport_hist_rcp45_branch_1500[:156], '-', color = 'k', linewidth = 0.5)

ax.set_xlim(0, 1200)
ax.set_ylim(-2, 22)
ax.set_xlabel('Model year after branching point')
ax.set_ylabel('Volume transport (Sv)')
ax.grid()

graphs	      	= graph_1 + graph_2 + graph_3 + graph_4 + graph_5
legend_labels 	= [l.get_label() for l in graphs]
legend_1	= ax.legend(graphs, legend_labels, loc='upper left', ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_title('b) AMOC strength at 26$^{\circ}$N ($F_H = 0.45$ Sv)')

#-----------------------------------------------------------------------------------------

ax2 		= fig.add_axes([0.617, 0.62, 0.33, 0.20])

ax2.plot([0.18, 0.18], [-2, 22], '--k')
ax2.plot([0.45, 0.45], [-2, 22], '--k')
ax2.text(0.18, 23, '0.18', verticalalignment='bottom', horizontalalignment='center', color = 'k', fontsize=10)
ax2.text(0.45, 23, '0.45', verticalalignment='bottom', horizontalalignment='center', color = 'k', fontsize=10)

ax2.plot(time_QE * 0.0003, transport_QE, '-', color = 'r', linewidth = 0.25)
ax2.plot(0.66 * 2 - time_QE * 0.0003, transport_QE, '-', color = 'b', linewidth = 0.25)
ax2.set_xlim(0, 0.66)
ax2.set_ylim(-2, 35)
ax2.set_xlabel('$F_H$ (Sv)')
ax2.set_ylabel('AMOC (Sv)')
ax2.grid()
ax2.set_title('Quasi-equilibrium')

ax2.quiver(1000 * 0.0003, 18, 3, 0, scale = 30, color = 'r', zorder = 10, width = 0.009)
ax2.quiver(1200 * 0.0003, 5, -3, 0, scale = 30, color = 'b', zorder = 10, width = 0.009)

show()
