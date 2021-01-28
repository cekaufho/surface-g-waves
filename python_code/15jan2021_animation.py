#=========================================================================================================
# SURFACE WAVES EXPERIMENT
# By Johannes and Christine 
# Python Code for Visualization
# Jan 14, 2021
#=========================================================================================================

# %% Import Modules

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import xarray

# %% Set- Up

data = xarray.open_dataset("10jan2021_1.cdf")
depth = data["ht"]
P, T = data["surf_press"], data["temp"]
U, V, W = data["u"], data["v"], data["w"]

# %% Tank dimensions
x = range(0, 60)
y = range(0, 60)
z = range(0, 2)
null = np.zeros(60)

#%% Surface Pressure

fig = plt.figure(figsize=(6,3), dpi=100)
ax1 = fig.add_subplot(121)
ax1.set_xlabel('Surface Pressure')

def updateLine(step, ax1):
    ax1.cla()
    ax1.axvline(color="grey", linestyle="--")
    ax1.plot(P[step, :, 0], y, color="red", label="x = 01")
    ax1.plot(P[step, :, 59], y, color="blue", label="x = 60")
    ax1.set_xlabel('Surface Pressure')
    ax1.set_xlim(-1E-5, 1E-5)
    ax1.set_ylim(0, 60)
    ax1.legend(loc='upper right', fontsize=8)

#%% Animation

ax2 = fig.add_subplot(122)
ax2.set_xlabel('X-Direction, Zonal')
ax2.set_ylabel('Y-Direction, Meridional')

def updateContour(step, ax2):
    ax2.cla()
    ax2.pcolor(x, y, P[step], cmap="seismic", vmin=-1e-6, vmax=1e-6)
    ax2.set_xlabel('X-Direction, Zondal')
    ax2.set_ylabel('Y-Direction, Meridional')    

#%% Plot Together

def updateALL(step, ax1, ax2):
    a = updateLine(step, ax1)
    b = updateContour(step, ax2)
    return b, a

#animation = FuncAnimation(fig, updateALL, frames=50, fargs=(ax2, ax1, ))
#animation.save('surface_pressure.mp4', writer='ffmpeg', fps=4);