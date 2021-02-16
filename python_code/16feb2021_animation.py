#=========================================================================================================
# SURFACE WAVES EXPERIMENT
# By Johannes and Christine 
# Python Code for Visualization
#=========================================================================================================

import numpy as np
import matplotlib.pyplot as plt
import numpy.ma as ma
from scipy.ndimage.filters import gaussian_filter
from matplotlib.animation import FuncAnimation
import xarray

# %% Set- Up

data = xarray.open_dataset("16feb2021_23.cdf")
depth = data["ht"]
P = np.array(data.variables["surf_press"][:])

# %% Tank dimensions
np.set_printoptions(suppress=True)

# Length in [km]
X = np.array(data.variables["xu"][:])
Y = np.array(data.variables["yu"][:])
Z = np.array(data.variables["zu"][:])

# Convert to [m] - Z already in [m]
X = X*1000.
Y = Y*1000.

# %% Print Statements

print("======================================================================")
print("The length of X in meters is ", X[-1], "meters with ", len(X), "grid boxes")
print("The length of Y in meters is ", Y[-1], "meters with ", len(Y), "grid boxes")
print("The length of Z in meters is ", Z[0], "meters with ", len(Z), "grid boxes")
print("The simulation runs for ", len(P), "seconds")
print("======================================================================")

#%% Surface Pressure

fig = plt.figure(figsize=(6,3), dpi=100)
ax1 = fig.add_subplot(121)
ax1.set_xlabel(r"Surface Pressure [$m^{2}$/$s^{2}$]")

def updateLine(step, ax1):
    ax1.cla()
    ax1.axvline(color="grey", linestyle="--")
    ax1.plot(P[step, :, 1], Y, color="red", label="x = beginning")
    ax1.plot(P[step, :, -1]-np.mean(P[step, :, -1]), Y, color="blue", label="x = end")
    ax1.set_xlabel(r"Surface Pressure [$m^{2}$/$s^{2}$]")
    ax1.set_xlim(-1E-8, 1E-8)
    ax1.set_ylim(Y[0], Y[-1])
    ax1.legend(loc='upper right', fontsize=8)

#%% Animation

ax2 = fig.add_subplot(122)
ax2.set_xlabel('Zonal, [$m$]')
ax2.set_ylabel('Meridional [$m$]')

def updateContour(step, ax2):
    ax2.cla()
    P_smooth = gaussian_filter(P[step], sigma = 3)
    P_mask = ma.masked_where(P[step] == 0, P_smooth)
    ax2.pcolor(X, Y, P_mask, vmin = -1E-8, vmax = 1E-8, cmap='seismic')
    ax2.set_xlabel('Zonal, [$m$]')
    ax2.set_ylabel('Meridional [$m$]')

#%% Plot Together

def updateALL(step, ax1, ax2):
    a = updateLine(step, ax1)
    b = updateContour(step, ax2)
    return b, a

animation = FuncAnimation(fig, updateALL, frames=len(P), fargs=(ax2, ax1, ))
animation.save('surface_pressure.mp4', writer='ffmpeg', fps=3);