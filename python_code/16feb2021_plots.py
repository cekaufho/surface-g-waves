#=========================================================================================================
# SURFACE WAVES EXPERIMENT
# By Johannes and Christine 
# Python Code for Visualization
#=========================================================================================================

import numpy as np
import matplotlib.pyplot as plt
import numpy.ma as ma
from scipy.ndimage.filters import gaussian_filter
import xarray

# %% Set- Up

data = xarray.open_dataset("16feb2021_26.cdf")
depth = data["ht"]
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
print("The simulation runs for ", len(P), "seconds")
print("======================================================================")

# %% Plot

for i in range(0, len(P), 1):    

    fig = plt.figure(figsize=(8,4))    
    #P_smooth = gaussian_filter(P[i], sigma=3)
    #P_mask = ma.masked_where(P[i] == 0, P_smooth)
    
    cs = plt.pcolor(X, Y, P[i], vmin = -1E-6, vmax = 1E-6, cmap='seismic')
    cbar = plt.colorbar(cs, label=r"Surface Pressure [$m^{2}$/$s^{2}$]") 
    
    plt.xlabel("Zonal [m]", fontsize=10)
    plt.ylabel("Meridional [m]", fontsize=10)
    plt.title("Surface Pressure vs. Time = " + str(i) + "s", fontsize=12)
    plt.show()