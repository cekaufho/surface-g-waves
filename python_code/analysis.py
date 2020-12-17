import numpy as np
import matplotlib.pylab as plt
import xarray
# import os

# Programmed by Christine Kaufhold, IfM 2020

# %% Set- Up

data = xarray.open_dataset("basin_3d_pyOM.cdf")
temp = data["temp"] # Temperature(t,z,y,x)
U = data["u"]

# %% Tank dimensions
x = np.arange(60)*0.01
y = np.arange(60)*0.01
z = np.arange(20)*0.01

# %% Plot T(x,y) at z = 10cm
    
for t in range(5):
    plt.figure()
    print("Time Step " + str(t+1) + "/100")
    plt.pcolormesh(x,y,U[t,10,:,:])
    #plt.colorbar(label='T [°C]')
    plt.title("temperature in x-y plane at z = 3cm")
    plt.xlabel("x [m]")
    plt.ylabel("y [m]")
    #plt.savefig('img/' + file_name, format='png')
    plt.show()
    plt.close
