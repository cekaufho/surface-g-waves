#=========================================================================================================
# SURFACE WAVES EXPERIMENT
# By Johannes and Christine 
# Python Code for Visualization
# December 17, 2020
#=========================================================================================================

# %% Import Modules

import numpy as np
import matplotlib.pylab as plt
import matplotlib.ticker as ticker
import seaborn as sns
import xarray
# import os

# %% Set- Up

data = xarray.open_dataset("17dec2020_1.cdf")
temp = data["temp"] # Temperature(t,z,y,x)
depth = data["ht"]
p = data["surf_press"]
U = data["u"]

# %% Tank dimensions
#x = np.arange(60)*0.1
#y = np.arange(60)*0.1
#z = np.arange(20)*0.1

# %% Plot T(x,y) at z = 10cm

plt.figure(figsize=(10,5))       
ax = sns.heatmap(p[-1], annot=False, fmt=".1f", cmap="viridis_r", vmin=0, vmax=0.2, cbar_kws={'label': 'Colourbar Kwargs'}) # annot = True gives actual depth numbers for the boxes

# Tick Frequency
ax.xaxis.set_major_locator(ticker.MultipleLocator(5))
ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
ax.yaxis.set_major_locator(ticker.MultipleLocator(10))
ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
ax.invert_yaxis()

plt.xlabel("X (Zonal) [cm]", fontsize=13)
plt.ylabel("Y (Meridional) [cm]", fontsize=13)
plt.title("Surface Gravity Waves - Double Slit", fontsize=15)
plt.show()