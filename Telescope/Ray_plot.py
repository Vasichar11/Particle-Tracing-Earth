import matplotlib.pyplot as plt
import numpy as np
import h5py

filename = "interpolated_ray_pwr1.000000.h5"

f1 = h5py.File("h5files/" + filename,"r")
lat = f1["lat_int"][()]
kx_ray = f1["kx_ray"]                    
kz_ray = f1["kz_ray"]
kappa_ray = f1["kappa_ray"]
Ezw   = f1["Ezw"]  
Bw_ray = f1["Bw_ray"]
w1 = f1["w1"]
w2 = f1["w2"]
R1 = f1["R1"]
R2 = f1["R2"]
Bzw = f1["Bzw"]

fig, ax = plt.subplots()
ax.plot(Bzw,)
ax.plot(Ezw,)
fig.savefig("simulationMM/Ray_plot" + filename + ".png")