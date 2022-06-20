import matplotlib.pyplot as plt
import numpy as np
import h5py

f1 = h5py.File("h5files/interpolated_ray_pwr10.h5","r")
lat = f1["lat_int"][()]




fig, ax = plt.subplots()

ax.plot()