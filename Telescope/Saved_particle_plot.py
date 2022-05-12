
import matplotlib.pyplot as plt
import h5py
import numpy as np

f2 = h5py.File("h5files/1000p_both.h5","r")
saved_id      = f2["saved_id"][()]
saved_deta_dt   = f2["saved_deta_dt"][()]
saved_alpha   = f2["saved_alpha"][()]
saved_time    = f2["saved_time"][()] 
f2.close()

#Sort to lists that reference each other
saved_alpha, saved_deta_dt = zip(*sorted(zip(saved_alpha, saved_deta_dt)))
fix,ax = plt.subplots()
ax.plot(np.rad2deg(saved_alpha),saved_deta_dt)
ax.set(xlabel="P.A",ylabel="deta/dt",title="alpha-deta/dt plot")
plt.savefig("alpha_detadt.png")