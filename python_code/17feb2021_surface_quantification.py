#=========================================================================================================
# SURFACE WAVES EXPERIMENT
# By Johannes and Christine 
# Python Code for Visualizationn
#=========================================================================================================

import numpy as np
import matplotlib.pyplot as plt
import numpy.ma as ma
from scipy import interpolate
from scipy.signal import savgol_filter
from scipy.ndimage.filters import gaussian_filter
import xarray

# %% Set- Up

data = xarray.open_dataset("13jan2021_1.cdf")
#data = xarray.open_dataset("17feb2021_03.cdf")
depth, time = data["ht"], np.array(data.variables['Time'][:])
P = np.array(data.variables["surf_press"][:])

# %% Tank dimensions
#np.set_printoptions(suppress=True)

# Length in [km]
X = np.array(data.variables["xu"][:])
Y = np.array(data.variables["yu"][:])

# Convert to [m] - Z already in [m]
X = X*1000.
Y = Y*1000.

# Length in [m]
Z = np.array(data.variables["zu"][:])

# %% Print Statements

print("======================================================================")
print("The length of X in meters is ", X[-1], "meters with ", len(X), "grid boxes")
print("The length of Y in meters is ", Y[-1], "meters with ", len(Y), "grid boxes")
print("The length of Z in meters is ", Z[0], "meters with ", len(Z), "grid boxes")
print("The simulation runs for ", len(P), "steps")
print("======================================================================")

# %% Take mean over all time

fig = plt.figure(figsize=(8,4))    
pressure = []

for i in range(9, 10, 1):
    P_i = P[i, :, -1] - np.mean(P[i, :, -1])
    pressure.append(P_i)

pressure = np.array(pressure)
pressure[(pressure>-1E-11) & (pressure<1E-11)] = 0

pressure_mean = np.mean(pressure, axis=0)
pressure_filter = gaussian_filter(pressure_mean, sigma=4)

# %% Plot

plt.plot(pressure_filter, Y, "k--", label=r"Gaussian Filter, $\sigma=4$", linewidth=1)
plt.scatter(pressure_mean, Y, s=15, color='c')
plt.vlines(0, Y[0], Y[-1], linewidth=1, color="grey", linestyle="--")
plt.xlabel(r"Surface Pressure [$m^{2}$/$s^{2}$]", fontsize=10)
plt.ylabel(r"Meridional [$m$]", fontsize=10)
plt.title("Mean Surface Pressure", fontsize=12)
plt.legend(loc="upper right")    
plt.show()