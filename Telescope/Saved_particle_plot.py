
import matplotlib.pyplot as plt
import h5py
import numpy as np

f2 = h5py.File("h5files/1000p_both.h5","r")
saved_id        = f2["saved_id"][()]
saved_deta_dt   = f2["saved_deta_dt"][()]
saved_lamda     = f2["saved_lamda"][()]
saved_Ekin      = f2["saved_Ekin"][()]
f2.close()

print(saved_id)
resonant = []
for i,j in zip(saved_deta_dt,saved_id):
    if abs(i)<0.1:
        resonant.append(j)


print(resonant)

#For every particle in the population the state for minimum deta_dt is saved. 
fix,ax = plt.subplots()
cc = ax.scatter(np.rad2deg(saved_lamda),saved_deta_dt,s=1, c=saved_Ekin,cmap='jet')
ax.set(xlabel="latitude[deg]",ylabel="deta_dt")
cbar=plt.colorbar(cc)
cbar.ax.set_title('Ekin[keV]')
plt.savefig("latitude-deta_dt-Ekin.png")

