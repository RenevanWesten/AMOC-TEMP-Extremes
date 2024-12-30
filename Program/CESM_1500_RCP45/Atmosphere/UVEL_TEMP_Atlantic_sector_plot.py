#Program plots the wind speeds and temperature over the Atlantic sector (45W - 15E)

from pylab import *
import numpy
import time
import glob, os
import math
import netCDF4 as netcdf
import matplotlib.colors as colors
from cartopy import crs as ccrs, feature as cfeature
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy import stats

#Making pathway to folder with all data
directory_1500_PI	= '../../../Data/CESM_1500_PI/'
directory_1500_RCP	= '../../../Data/CESM_1500_RCP45/'

def Welch(data_1, data_2):
    """Conducts the Welch t-test"""
	
    #Determine the means
    mean_1	= np.mean(data_1)
    mean_2	= np.mean(data_2)
	
    #Determine the corrected sample standard deviations
    std_1	= np.sqrt(1.0 / (len(data_1) - 1) * np.sum((data_1 - mean_1)**2.0))
    std_2	= np.sqrt(1.0 / (len(data_2) - 1) * np.sum((data_2 - mean_2)**2.0))

    #Determine the Welch t-value
    t_welch	= (mean_1 - mean_2) / np.sqrt((std_1**2.0 / len(data_1)) + (std_2**2.0 / len(data_2)))

    #Determine the degrees of freedome (dof)
    dof	= ((std_1**2.0 / len(data_1)) + (std_2**2.0 / len(data_2)))**2.0 / ((std_1**4.0 / (len(data_1)**2.0 * (len(data_1) - 1))) + (std_2**4.0 / (len(data_2)**2.0 * (len(data_2) - 1))))

    #Get the significance levels and the corresponding critical values (two-sided)
    sig_levels 	= np.arange(50, 100, 0.5) / 100.0
    t_crit		= stats.t.ppf((1.0 + sig_levels) / 2.0, dof)

    #Get the indices where the significance is exceeding the critical values
    sig_index	= np.where(fabs(t_welch) > t_crit)[0]
    significant	= 0.0

    if len(sig_index) > 0:
        #If there are significance values, take the highest significant level
        significant = sig_levels[sig_index[-1]]

    return significant

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

month_start	= 1
month_end	= 1

#-----------------------------------------------------------------------------------------

fh 		    = netcdf.Dataset(directory_1500_RCP+'Atmosphere/UVEL_TEMP_Atlantic_sector_month_'+str(month_start)+'-'+str(month_end)+'.nc', 'r')

pres		= fh.variables['lev'][:]
lat		    = fh.variables['lat'][:]				
u_vel		= fh.variables['U'][:]			#Zonal velocity (m / s)	
temp		= fh.variables['TEMP'][:]		#Potential temperature (K)

fh.close()

#-----------------------------------------------------------------------------------------

fh          = netcdf.Dataset(directory_1500_PI+'Atmosphere/UVEL_TEMP_Atlantic_sector_month_'+str(month_start)+'-'+str(month_end)+'.nc', 'r')
		
u_vel_ref	= fh.variables['U'][:]			#Zonal velocity (m / s)	
temp_ref	= fh.variables['TEMP'][:]		#Potential temperature (K)

fh.close()

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

#Determine the position of het subpolar jet (below 100 hPa and above 40N)
pres_min_index	= (np.abs(pres - 140)).argmin()
lat_min_index	= (np.abs(lat - 40)).argmin()

lat_jet		= ma.masked_all(len(u_vel))
pres_jet	= ma.masked_all(len(u_vel))
lat_jet_ref	= ma.masked_all(len(u_vel))
pres_jet_ref	= ma.masked_all(len(u_vel))

for year_i in range(len(u_vel)):
    #Get the maximum of the subpolar jet
    jet_max_index	= np.where(u_vel[year_i, pres_min_index:, lat_min_index:] == np.max(u_vel[year_i, pres_min_index:, lat_min_index:]))
    jet_pres_index	= jet_max_index[0][0] + pres_min_index
    jet_lat_index	= jet_max_index[1][0] + lat_min_index
	
    lat_jet[year_i]	= lat[jet_lat_index]
    pres_jet[year_i]= pres[jet_pres_index]

    #Get the maximum of the subpolar jet (for the reference)
    jet_max_index	= np.where(u_vel_ref[year_i, pres_min_index:, lat_min_index:] == np.max(u_vel_ref[year_i, pres_min_index:, lat_min_index:]))
    jet_pres_index	= jet_max_index[0][0] + pres_min_index
    jet_lat_index	= jet_max_index[1][0] + lat_min_index
	
    lat_jet_ref[year_i]	= lat[jet_lat_index]
    pres_jet_ref[year_i]	= pres[jet_pres_index]

print('Jet latitudinal position:', np.mean(lat_jet))
print('Jet latitudinal position (ref):', np.mean(lat_jet_ref))
print('Welch t-test for significance:', Welch(lat_jet, lat_jet_ref))

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

u_vel_plot	= np.mean(u_vel, axis = 0) - np.mean(u_vel_ref, axis = 0)

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

CS	= ax.contourf(lat, pres, u_vel_plot, levels = np.arange(-10, 10.1, 1), extend = 'both', cmap = 'PuOr_r')
cbar	= colorbar(CS, ticks = np.arange(-10, 10.1, 5))
cbar.set_label('Zonal velocity difference (m s$^{-1}$)')

CS2	= ax.contour(lat, pres, -np.mean(u_vel_ref, axis = 0), levels = [-30], colors = 'gray', linewidths = 2, zorder = 10)
CS2	= ax.contour(lat, pres, np.mean(temp, axis = 0), levels = [300], colors = 'firebrick', linewidths = 2, zorder = 10)
CS2	= ax.contour(lat, pres, np.mean(temp, axis = 0), levels = [320], colors = 'royalblue', linewidths = 2, zorder = 10)
CS2	= ax.contour(lat, pres, np.mean(temp, axis = 0), levels = [350], colors = 'k', linewidths = 2, zorder = 10)

ax.set_ylabel('Pressure (hPa)')
ax.set_xlim(0, 80)
ax.set_ylim(950, 50)
ax.set_yscale('log')

ax.set_xticks(np.arange(0, 80.1, 20))
ax.set_xticklabels(['Eq', '20$^{\circ}$N', '40$^{\circ}$N', '60$^{\circ}$N', '80$^{\circ}$N'])
ax.set_yticks([50, 100, 200, 300, 400, 500, 600, 850])
ax.set_yticklabels([50, 100, 200, 300, '', 500, '', 850])

graph_1		= ax.plot([90, 90], [-1, -1], linestyle = '-', color = 'k', linewidth = 2, label = r'$\theta$ = 350 K')
graph_2		= ax.plot([90, 90], [-1, -1], linestyle = '-', color = 'royalblue', linewidth = 2, label = r'$\theta$ = 320 K')
graph_3		= ax.plot([90, 90], [-1, -1], linestyle = '-', color = 'firebrick', linewidth = 2, label = r'$\theta$ = 300 K')
graph_4		= ax.plot([90, 90], [-1, -1], linestyle = '--', color = 'gray', linewidth = 2, label = '30 m s$^{-1}$ (AMOC on)')

graphs		= graph_1 + graph_2 + graph_3 + graph_4
legend_labels = [l.get_label() for l in graphs]
legend_1      = ax.legend(graphs, legend_labels, loc = 'upper left', ncol=1, numpoints = 1, framealpha = 1.0).set_zorder(20)

ax.set_title('l) Zonal velocity, RCP4.5 ($F_H = 0.45$ Sv)')

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

lat_index_1	= (np.abs(lat - 45)).argmin()
lat_index_2	= (np.abs(lat - 55)).argmin()+1
lat_index_3	= (np.abs(lat - 50)).argmin()

fig, ax	= subplots()

CS	= ax.contourf(lat, pres, np.mean(temp, axis = 0) - np.mean(temp_ref, axis = 0), levels = np.arange(-10, 10.1, 1), extend = 'both', cmap = 'RdBu_r')
cbar	= colorbar(CS, ticks = np.arange(-10, 10.1, 5))
cbar.set_label('Potential temperature difference (K)')

CS2	= ax.contour(lat, pres, np.mean(temp, axis = 0), levels = [300], colors = 'firebrick', linewidths = 2, zorder = 10)
CS2	= ax.contour(lat, pres, np.mean(temp, axis = 0), levels = [320], colors = 'royalblue', linewidths = 2, zorder = 10)
CS2	= ax.contour(lat, pres, np.mean(temp, axis = 0), levels = [350], colors = 'k', linewidths = 2, zorder = 10)

ax.set_ylabel('Pressure (hPa)')
ax.set_xlim(0, 70)
ax.set_ylim(950, 50)
ax.set_yscale('log')

ax.set_xticks(np.arange(0, 80.1, 20))
ax.set_xticklabels(['Eq', '20$^{\circ}$N', '40$^{\circ}$N', '60$^{\circ}$N', '80$^{\circ}$N'])
ax.set_yticks([50, 100, 200, 300, 400, 500, 600, 850])
ax.set_yticklabels([50, 100, 200, 300, '', 500, '', 850])

graph_1		= ax.plot([90, 90], [-1, -1], linestyle = '-', color = 'k', linewidth = 2, label = r'$\theta$ = 350 K')
graph_2		= ax.plot([90, 90], [-1, -1], linestyle = '-', color = 'royalblue', linewidth = 2, label = r'$\theta$ = 320 K')
graph_3		= ax.plot([90, 90], [-1, -1], linestyle = '-', color = 'firebrick', linewidth = 2, label = r'$\theta$ = 300 K')

graphs		= graph_1 + graph_2 + graph_3
legend_labels = [l.get_label() for l in graphs]
legend_1      = ax.legend(graphs, legend_labels, loc = 'upper left', ncol=1, numpoints = 1, framealpha = 1.0).set_zorder(20)

fig.subplots_adjust(top=0.850)

fig.suptitle('k) Potential temperature, RCP4.5 ($F_H = 0.45$ Sv)')

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

temp_grad	= ma.masked_all(len(pres))
TWB		    = ma.masked_all(len(pres))
    
P_lev_diff	= np.asarray([3, 3, 4 , 5, 7.5, 10, 10, 10, 15, 20, 20,   20,  20,  20,  20,  35, 50 ,  50,  50,  75, 100, 100, 100, 75 , 50,37.5, 25,   25, 25, 25])

factor_1	= np.log((1.0 + np.e)/2)
factor_2	= 1.0 - factor_1

for pres_i in range(len(pres)):
    #Loop over each pressure level and determine the trend
    grid_y	= lat[lat_index_1:lat_index_2]
    grid_y	= grid_y - grid_y[0]
	
    #Convert to meter
    grid_y	= grid_y * 6371000 * 2 * np.pi / 360.0

    a, b	= np.polyfit(grid_y, np.mean(temp[:, pres_i, lat_index_1:lat_index_2], axis = 0), 1)
	
    #Save the meridional temperature gradient
    temp_grad[pres_i]	= a
	
    #Determine the thermal wind balance
    if pres_i == 0:
        TWB[pres_i]	= 0.
		
    else:
        u_vel_diff_1	= 8.314 / (2 * 7.29 * 10**(-5.0) * np.sin(lat[lat_index_3] * 2 * np.pi / 360.)) * temp_grad[pres_i-1]
        u_vel_diff_2	= 8.314 / (2 * 7.29 * 10**(-5.0) * np.sin(lat[lat_index_3] * 2 * np.pi / 360.)) * temp_grad[pres_i]
		
        TWB[pres_i]	= TWB[pres_i-1] + u_vel_diff_1 * np.log(P_lev_diff[pres_i-1]) * factor_2 + u_vel_diff_2 * np.log(P_lev_diff[pres_i]) * factor_1

#Take the 850-hPa level as a reference
TWB	= TWB - TWB[24] + np.mean(u_vel[:, 24, lat_index_3], axis = 0)

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

graph_1	= ax.plot(np.mean(u_vel[:, :, lat_index_3], axis = 0), pres, '-', color = 'k', linewidth = 2.0, label = '$u$ (50$^{\circ}$N)')
graph_2	= ax.plot(TWB, pres, '-', color = 'royalblue', linewidth = 2.0, label = 'TWB')

ax.set_xlabel('Zonal velocity (m s$^{-1}$)')
ax.set_ylabel('Pressure (hPa)')
ax.set_xlim(0, 40)
ax.set_ylim(950, 50)
ax.set_yscale('log')
ax.grid()

ax2 	 = ax.twiny()
ax2.set_xlim(0, 50)

graph_3	= ax2.plot(temp_grad * 1000000, pres, '-', linewidth = 2.0, color = 'firebrick', label = '$\partial_y T$ (45$^{\circ}$N - 55$^{\circ}$N)')
ax2.set_xlim(-20, 20)
ax2.set_xlabel('Meridional temperature gradient (K per 1000 km)')

ax.set_yticks([50, 100, 200, 300, 400, 500, 600, 850])
ax.set_yticklabels([50, 100, 200, 300, '', 500, '', 850])

graphs		= graph_1 + graph_2 + graph_3
legend_labels = [l.get_label() for l in graphs]
legend_1      = ax.legend(graphs, legend_labels, loc = 'upper left', ncol=1, numpoints = 1, framealpha = 1.0).set_zorder(20)

fig.subplots_adjust(top=0.850)

fig.suptitle('l) Thermal wind balance, RCP4.5 ($F_H = 0.45$ Sv)')

show()

