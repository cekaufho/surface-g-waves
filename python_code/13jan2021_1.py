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

"""
for i in range(0, 2000, 1):    

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
    
# Tick Frequency
#ax.xaxis.set_major_locator(ticker.MultipleLocator(5))
#ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
#ax.yaxis.set_major_locator(ticker.MultipleLocator(10))
#ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
"""


#%% test

fig = plt.figure()
ax = fig.add_subplot(111)
plt.plot(y, P[0, :, 59])
plt.ylim(-1E-5, 1E-5)
plt.xlabel('Y-Direction, Meridional', fontsize=12)
plt.ylabel('Surface Pressure', fontsize=12)
plt.title("Surface Pressure vs. Y-Direction", fontsize=13)

def create_frame(step, ax):
    ax.cla()
    plt.plot(y, P[step, :, 59])
    plt.ylim(-1E-5, 1E-5)
    plt.xlabel('Y-Direction, Meridional', fontsize=12)
    plt.ylabel('Surface Pressure', fontsize=12)
    plt.title("Surface Pressure vs. Y-Direction", fontsize=13)

from matplotlib.animation import FuncAnimation
animation = FuncAnimation(fig, create_frame, frames=100, fargs=(ax,))
animation.save('animation_test.mp4', writer='ffmpeg', fps=5);


#%% Animate

#fig = plt.figure(figsize=(8,4))
#ax = fig.gca()
#surf_press = sns.heatmap((P[1]), ax=ax, annot=False, fmt=".1f", cmap="viridis")# vmin=-1e-7, vmax=1e-7)
#ax.set_title("Surface Pressure")


"""
fig = plt.figure()
ax = fig.add_subplot(111)
plt.pcolor(x, y, P[0], cmap="seismic", vmin=-1e-6, vmax=1e-6)
plt.colorbar()
plt.xlabel('Zondal', fontsize=12)
plt.ylabel('Meridional', fontsize=12)
plt.title("Surface Pressure", fontsize=13)

def create_frame(step, ax):
    ax.cla()
    plt.pcolor(x, y, P[step], cmap="seismic", vmin=-1e-6, vmax=1e-6)
    #sns.heatmap((P[step]), ax=ax, annot=False, fmt=".1f", cmap="viridis")
    #sns.lineplot(x[:step], y[:step], ax=ax)

from matplotlib.animation import FuncAnimation
#fig = plt.figure()
#ax = fig.gca()
#create_frame(10, ax)
animation = FuncAnimation(fig, create_frame, frames=500, fargs=(ax,))
animation.save('animation_test.mp4', writer='ffmpeg', fps=5);

#import os
#os.system("ffmpeg -i C:\\Users\\owner\\Desktop\\Processes Code\\animation_test.mp4 C:\\Users\\owner\\Desktop\\Processes Code\\animation_test.gif")

#animation = FuncAnimation(fig, create_frame, fargs=(ax,), frames=10, interval=200, blit=True)
#animation.save('0001_11jan2021.mp4', fps=20)
"""