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
directory	= '../../../Program/'
	
#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

experiment     = ['CESM_0600_PI', 'CESM_3800_PI', 'CESM_0600_RCP45', 'CESM_0600_RCP85', 'CESM_1500_PI', 'CESM_2900_PI', 'CESM_1500_RCP45', 'CESM_1500_RCP85']

for exp_i in experiment:
    #Loop over each experiment
    os.system('python '+directory+exp_i+'/Atmosphere/Jet_200_hPa_Atlantic_sector_plot.py')