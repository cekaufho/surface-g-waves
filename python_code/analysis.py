#=========================================================================================================
# SURFACE WAVES EXPERIMENT
# By Johannes and Christine 
# Python Code for Visualization
#=========================================================================================================

# %% Import Modules

import numpy as np
import matplotlib.pylab as plt
import seaborn as sns
import xarray
# import os

# %% Set- Up

data = xarray.open_dataset("2slit.cdf")
temp = data["temp"] # Temperature(t,z,y,x)
depth = data["ht"]
p=data["surf_press"]
U = data["u"]

# %% Tank dimensions
x = np.arange(60)*0.1
y = np.arange(60)*0.1
z = np.arange(20)*0.1

# %% Plot T("x,y) at z = "10cm
f, ax = plt.subplots(figsize=(18, 12))
sns.heatmap(depth,annot=True,fmt=".1f", ax=ax,cmap="mako_r")


"""
# %% Plot T(x,y) at z = 10cm
    
for t in range(5):
    plt.figure()
    print("Time Step " + str(t+1) + "/100")
    plt.pcolormesh(x,y,U[t,10,:,:])
    #plt.colorbar(label='T [°C]')
    plt.pcolormesh(depth,cmap="viridis", rasterized=True)
    #plt.savefig('img/' + file_name, format='png')
    plt.title("temperature in x-y plane at z = 3cm")
    plt.xlabel("x [m]")
    plt.ylabel("y [m]")
    #plt.savefig('img/' + file_name, format='png')
    plt.show()
    plt.close
"""