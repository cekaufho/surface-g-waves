#=========================================================================================================
# SURFACE WAVES EXPERIMENT
# By Johannes and Christine 
# Python Code for Visualization
# Jan 14, 2021
#=========================================================================================================

# %% Import Modules

import numpy as np
#import matplotlib.pylab as plt
import matplotlib.pyplot as plt
#import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
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

#%% Surface Pressure

fig = plt.figure(figsize=(6,3), dpi=100)
ax1 = fig.add_subplot(121)
#plt.plot(y, P[0, :, 59])
#plt.ylim(-1E-5, 1E-5)
#plt.xlabel('Y-Direction, Meridional', fontsize=12)
#plt.ylabel('Surface Pressure', fontsize=12)
#plt.title("Surface Pressure", fontsize=13)
ax1.set_xlabel('Surface Pressure')
ax1.set_ylabel('Y-Direction')

def updateLine(step, ax1):
    ax1.cla()
    ax1.plot(P[step, :, 59], y)
    ax1.set_xlabel('Surface Pressure')
    ax1.set_xlim(-1E-5, 1E-5)
    ax1.set_ylim(0, 60)
    #plt.ylim(-1E-5, 1E-5)
    #plt.xlabel('Y-Direction, Meridional', fontsize=12)
    #plt.ylabel('Surface Pressure', fontsize=12)
    #plt.title("Surface Pressure vs. Y-Direction", fontsize=13)

#animation = FuncAnimation(fig, create_frame, frames=10, fargs=(ax,))
#animation.save('surf_press_y.mp4', writer='ffmpeg', fps=5);

#%% Animation

ax2 = fig.add_subplot(122)
#plt.pcolor(x, y, P[0], cmap="seismic", vmin=-1e-6, vmax=1e-6)
#plt.colorbar()
#plt.xlabel('Zondal', fontsize=12)
#plt.ylabel('Meridional', fontsize=12)
#plt.title("Surface Pressure", fontsize=13)
ax2.set_xlabel('X-Direction')
ax2.set_ylabel('Y-Direction')

def updateContour(step, ax2):
    ax2.cla()
    ax2.pcolor(x, y, P[step], cmap="seismic", vmin=-1e-6, vmax=1e-6)
    ax2.set_xlabel('X-Direction')
    ax2.set_ylabel('Y-Direction')    

#animation = FuncAnimation(fig, create_frame, frames=10, fargs=(ax,))
#animation.save('surf_press_xy.mp4', writer='ffmpeg', fps=5);

#%% Plot Together

def updateALL(step, ax1, ax2):
    a = updateLine(step, ax1)
    b = updateContour(step, ax2)
    return b, a

animation = FuncAnimation(fig, updateALL, frames=1000, fargs=(ax2, ax1, ))
animation.save('surface_pressure.mp4', writer='ffmpeg', fps=5);