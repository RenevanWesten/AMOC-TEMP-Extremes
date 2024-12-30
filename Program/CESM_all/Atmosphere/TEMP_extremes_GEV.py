#Program determines the local GEV fits for both the minima and maxima temperatures

from pylab import *
import numpy
import time
import glob, os
import math
import netCDF4 as netcdf
import matplotlib.colors as colors
from scipy import stats
from scipy.stats import genextreme
from matplotlib.colors import LogNorm
from scipy.stats import chisquare

#Making pathway to folder with all data
directory	= '../../../Data/'

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

experiment     = ['CESM_0600_PI', 'CESM_3800_PI', 'CESM_0600_RCP45', 'CESM_0600_RCP85', 'CESM_1500_PI', 'CESM_2900_PI', 'CESM_1500_RCP45', 'CESM_1500_RCP85']

for exp_i in experiment:
    #Loop over each experiment
    
    for temp_extreme_type in ['minima', 'maxima']:
        #Loop over the two extreme types

        #-----------------------------------------------------------------------------------------
        
        fh 		= netcdf.Dataset(directory+exp_i+'/Atmosphere/TEMP_2m_extremes.nc', 'r')
        
        time		= fh.variables['time'][:] 	
        lat		    = fh.variables['lat'][:] 			
        lon		    = fh.variables['lon'][:] 			
        temp_min	= fh.variables['TEMP_year_min'][:] 
        temp_max	= fh.variables['TEMP_year_max'][:] 
        
        fh.close()

        #-----------------------------------------------------------------------------------------
        
        shape_all	= ma.masked_all((len(lat), len(lon)))
        loc_all		= ma.masked_all((len(lat), len(lon)))
        scale_all	= ma.masked_all((len(lat), len(lon)))
        p_value_all	= ma.masked_all((len(lat), len(lon)))
        
        #-----------------------------------------------------------------------------------------

        for lat_i in range(len(lat)):
            print('Experiment:', exp_i)
            print('Extreme type:', temp_extreme_type)
            print('Latitude:', lat[lat_i])
            print()
            for lon_i in range(len(lon)):

                if temp_extreme_type == 'maxima':
                    #Get the data for each grid cell for the yearly maxima
                    data	= temp_max[:, lat_i, lon_i]	

                if temp_extreme_type == 'minima':
                    #Get the data for each grid cell for the yearly minima, add minus sign
                    data	= -temp_min[:, lat_i, lon_i]	

                #Fit GEV distribution to stationary data (give along the initial guess, otherwise wrong fits)
                shape_1, loc_1, scale_1	= genextreme.fit(data, loc = np.mean(data), scale = np.std(data))
                shape_2, loc_2, scale_2	= genextreme.fit(data, scale = np.std(data))
                shape_3, loc_3, scale_3	= genextreme.fit(data, loc = np.mean(data))
                shape_4, loc_4, scale_4	= genextreme.fit(data)

                #Reverse shape (Python reverses the shape for plotting)
                shape_1, shape_2, shape_3, shape_4	= -shape_1, -shape_2, -shape_3, -shape_4

                #Same binsize everywhere as the temperature range is sufficient
                bin_size = 0.2

                #Get the minimum and maximum range
                data_max	= int(np.max(data) / bin_size + 1) * bin_size

                #For the minimum range, select the 2.5% percentile
                #This percentile is used to exclude strong tails in the fit, which then result in a poor KS fit
                #While the fit is actually very accurate
                per_value	= 2.5
                data_min	= np.percentile(data, per_value)
                data_min	= int(data_min / bin_size) * bin_size

                bins		= np.arange(data_min, data_max + bin_size / 10, bin_size)
                data_bins	= np.zeros(len(bins))
                cdist_GEV_1	= genextreme.cdf(bins, -shape_1, loc_1, scale_1)
                cdist_GEV_2	= genextreme.cdf(bins, -shape_2, loc_2, scale_2)
                cdist_GEV_3	= genextreme.cdf(bins, -shape_3, loc_3, scale_3)
                cdist_GEV_4	= genextreme.cdf(bins, -shape_4, loc_4, scale_4)

                for data_i in range(len(data)):
                    #Get the histrogram of the data
                    index			= np.where(data[data_i] >= bins)[0]

                    if len(index) > 0:
                        index	= index[-1]
                    else:
                        continue

                    data_bins[index]	+= 1.0 / len(data)

                #Determine the cummulative distribution
                for bin_i in range(1, len(bins)):
                    data_bins[bin_i]	= data_bins[bin_i] + data_bins[bin_i-1]

                #The cummulative distribution should add up to 1.0 (due to cut-off percentile)
                data_bins	+= 1.0 - data_bins[-1]

                #Determine the KS-test (cumulative distributions, below 0.05 the distributions are not the same)
                bins		= bins[:-1]
                cdist_GEV_1	= cdist_GEV_1[:-1]
                cdist_GEV_2	= cdist_GEV_2[:-1]
                cdist_GEV_3	= cdist_GEV_3[:-1]
                cdist_GEV_4	= cdist_GEV_4[:-1]
                data_bins	= data_bins[:-1]

                KS_1, p_value_1	= stats.kstest(data_bins, cdist_GEV_1)
                KS_2, p_value_2	= stats.kstest(data_bins, cdist_GEV_2)
                KS_3, p_value_3	= stats.kstest(data_bins, cdist_GEV_3)
                KS_4, p_value_4	= stats.kstest(data_bins, cdist_GEV_4)

                #Select the lowest KS value (which has the highest p-value)
                index_min	= np.argmin([KS_1, KS_2, KS_3, KS_4])
		
                if index_min == 0:
                    #Save the data for the first fit
                    shape_all[lat_i, lon_i]	= shape_1
                    loc_all[lat_i, lon_i]	= loc_1
                    scale_all[lat_i, lon_i]	= scale_1
                    p_value_all[lat_i, lon_i]   = p_value_1

                if index_min == 1:
                    #Save the data for the second fit
                    shape_all[lat_i, lon_i]	= shape_2
                    loc_all[lat_i, lon_i]	= loc_2
                    scale_all[lat_i, lon_i]	= scale_2
                    p_value_all[lat_i, lon_i]   = p_value_2

                if index_min == 2:
                    #Save the data for the third fit
                    shape_all[lat_i, lon_i]	= shape_3
                    loc_all[lat_i, lon_i]	= loc_3
                    scale_all[lat_i, lon_i]	= scale_3
                    p_value_all[lat_i, lon_i]   = p_value_3

                if index_min == 3:
                    #Save the data for the fourth fit
                    shape_all[lat_i, lon_i]	= shape_4
                    loc_all[lat_i, lon_i]	= loc_4
                    scale_all[lat_i, lon_i]	= scale_4
                    p_value_all[lat_i, lon_i]   = p_value_4

                if temp_extreme_type == 'minima':
                    #Add minus sign for the relative GEV variables
                    scale_all[lat_i, lon_i]	= -scale_all[lat_i, lon_i]
                    loc_all[lat_i, lon_i]	= -loc_all[lat_i, lon_i]

        fh = netcdf.Dataset(directory+exp_i+'/Atmosphere/TEMP_2m_extremes_GEV_fit_'+temp_extreme_type+'.nc', 'w')
        fh.createDimension('lat', len(lat))
        fh.createDimension('lon', len(lon))

        fh.createVariable('lat', float, ('lat'), zlib=True)
        fh.createVariable('lon', float, ('lon'), zlib=True)
        fh.createVariable('shape', float, ('lat', 'lon'), zlib=True)
        fh.createVariable('loc', float, ('lat', 'lon'), zlib=True)
        fh.createVariable('scale', float, ('lat', 'lon'), zlib=True)
        fh.createVariable('p_value', float, ('lat', 'lon'), zlib=True)
        
        fh.variables['lat'].long_name		= 'Array of latitudes'
        fh.variables['lon'].long_name		= 'Array of longitudes'
        fh.variables['shape'].long_name		= 'Shape parameter (GEV)'
        fh.variables['loc'].long_name		= 'Location parameter (GEV)'
        fh.variables['scale'].long_name		= 'Scale parameter (GEV)'
        fh.variables['p_value'].long_name	= 'P-value of fit'
        
        fh.variables['lat'].units 		= 'degrees N'
        fh.variables['lon'].units 		= 'degrees E'
        
        #Writing data to correct variable	
        fh.variables['lat'][:]     		= lat
        fh.variables['lon'][:]     		= lon
        fh.variables['shape'][:] 		= shape_all
        fh.variables['loc'][:] 			= loc_all
        fh.variables['scale'][:] 		= scale_all
        fh.variables['p_value'][:] 		= p_value_all
        
        fh.close()
