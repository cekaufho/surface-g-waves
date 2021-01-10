#=========================================================================================================
# SURFACE WAVES EXPERIMENT
# By Johannes and Christine 
# Python Code for Visualization
# December 17, 2020
#=========================================================================================================

# %% Import Modules

#import numpy as np
import matplotlib.pylab as plt
#import matplotlib.animation as animation
#import matplotlib.ticker as ticker
from matplotlib.colors import ListedColormap
import seaborn as sns
import xarray
# import os

# %% Set- Up

data = xarray.open_dataset("02jan2021_1.cdf")
depth = data["ht"]
P, T = data["surf_press"], data["temp"]
U, V, W = data["u"], data["v"], data["w"]

# %% Tank dimensions
#x = np.arange(60)*0.1
#y = np.arange(60)*0.1
#z = np.arange(20)*0.1

# %% Plot

for i in range(0, 2000):    
    fig = plt.figure(figsize=(8,4))
    ax = sns.heatmap((P[-i]-0.02314), annot=False, fmt=".1f", cmap="viridis") #vmin=-1e-5, vmax=1e-5)
    
    plt.xlabel("Zonal", fontsize=13)
    plt.ylabel("Meridional", fontsize=13)
    #ax.invert_yaxis()

    #plt.savefig('00' + str(i) + '.png')
    plt.title("Surface Pressure", fontsize=15)
    plt.show()
    
# Tick Frequency
#ax.xaxis.set_major_locator(ticker.MultipleLocator(5))
#ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
#ax.yaxis.set_major_locator(ticker.MultipleLocator(10))
#ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
