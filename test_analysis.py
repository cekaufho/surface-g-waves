import numpy as np
import matplotlib.pylab as pl
import os
import xarray

# Programmed by Christine Kaufhold, IfM 2020

# %% Set- Up

x= xarray.open_dataset("pyOM.cdf")