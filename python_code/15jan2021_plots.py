#=========================================================================================================
# SURFACE WAVES EXPERIMENT
# By Johannes and Christine 
# Python Code for Visualization
# Jan 13, 2021
#=========================================================================================================

# %% Import Modules

import numpy as np
#import matplotlib.pylab as plt
import matplotlib.pyplot as plt
#import matplotlib.animation as animation
#import matplotlib.ticker as ticker
#from matplotlib.colors import ListedColormap
from mpl_toolkits.mplot3d import Axes3D
from pylab import xlabel, ylabel, title
import seaborn as sns
import xarray
# import os

# %% Set- Up

data = xarray.open_dataset("10jan2021_1.cdf")
depth = data["ht"]
P, T = data["surf_press"], data["temp"]
U, V, W = data["u"], data["v"], data["w"]

# %% Tank dimensions
x = range(0, 60)
y = range(0, 60)
z = range(0, 2)

# %% Plot

for i in range(0, 5, 1):    

    fig = plt.figure(figsize=(8,4))
    #ax = fig.add_subplot(111, projection='3d')
    #ax.invert_yaxis()
    #ax.plot_surface(x, y, P[i], cmap="jet")
    #ax.set_zlim(-1e-7, 1e-7)
    #ax.view_init(65, 65)
    
    ax = sns.heatmap((P[-i]), annot=False, fmt=".1f", cmap="viridis")# vmin=-1e-7, vmax=1e-7)
    #-0.02314
    
    plt.xlabel("Zonal", fontsize=13)
    plt.ylabel("Meridional", fontsize=13)

    #ax.set_xlabel('X-axis', fontsize=14)
    #ax.set_ylabel('Y-axis', fontsize=14)
    #ax.set_zlabel('Z-axis', fontsize=14)
    #ax.plot(solution[0], solution[1], solution[2], color="k")

    #plt.savefig('00' + str(i) + '.png')
    plt.title("Surface Pressure", fontsize=15)
    plt.show()
    