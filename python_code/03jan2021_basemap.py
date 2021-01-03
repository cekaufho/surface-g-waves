from netCDF4 import Dataset, num2date
from scipy.interpolate import griddata
from matplotlib.colors import ListedColormap
import matplotlib.pylab as plt
from matplotlib import cm
import numpy.ma as ma
import pandas as pd
import numpy as np
import xarray
import csv
import os

from datetime import datetime
dateparse = lambda x: datetime.strptime(x, '%Y-%m-%d %H:%M:%S')
os.environ["PROJ_LIB"] = "C:\\Users\\owner\\anaconda3\\Library\\share"
from mpl_toolkits.basemap import Basemap

#data_r = xarray.open_dataset("SMOD_1979-2017_v04r00.nc")
#data = Dataset("SMOD_1979-2017_v04r00.nc", mode='r')
#data_range = data.variables["range"][:]
#smod = data.variables["SMOD"][:]
#latitude = data.variables["latitude"][:]
#longitude = data.variables["longitude"][:]
#data.close()
data = xarray.open_dataset("02jan2021_5.cdf")
P = data["surf_press"]
data.close()

plt.figure(figsize=(7,4), dpi=100) 
#m = Basemap(projection='npstere',boundinglat=55,lon_0=10,resolution='l')
m = Basemap(projection='gall',llcrnrlat=53.6,urcrnrlat=54.1,\
            llcrnrlon=7.1,urcrnrlon=7.9,resolution='h')
    
m.drawmapboundary(fill_color='aqua')
m.drawcoastlines()
m.fillcontinents(color='#cc9955',lake_color='#cc9955')
m.drawparallels(np.arange(10,90,1), labels=[True,False,False,False])
m.drawmeridians(np.arange(-180,180,1), labels=[False,False,False,True])
m.drawmapscale(7.25, 54.05, 6, 54, 25, barstyle='fancy')

longitude = np.linspace(7.2, 7.9, 60)
latitude = np.linspace(53.7, 54.2, 60)

lon, lat = np.meshgrid(longitude, latitude)
rot  = 0
shift = 0.5 #shift coords

lon   =  np.cos(np.deg2rad(rot))*lon + np.sin(np.deg2rad(rot))*lon
lat   = -np.sin(np.deg2rad(rot))*lat + np.cos(np.deg2rad(rot))*lat

xi, yi = m(lon, lat)

yi = yi + 3

# Define colourbar
#jet_big = cm.get_cmap('jet', 512)
#newcmp = ListedColormap(jet_big(np.linspace(0., 1, 256)))

#cs = m.contourf(xi, yi, P[-1], 500, cmap="jet")
cs = m.pcolor(xi, yi, P[29].transpose(), cmap="jet")
plt.title(r"Frisian Island Chain", fontsize=15)
cbar = m.colorbar(cs, location='right', label="Sea Surface Pressure")
#smod_max = smod[20].max()
#plt.clim(0, smod_max)
plt.show()