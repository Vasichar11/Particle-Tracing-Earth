import numpy as np 
import matplotlib.pyplot as plt
from matplotlib import cm
import h5py


D2R=np.pi/180
R2D=1/D2R


#Reading the data from the file and then closing it.
f = h5py.File("./h5files/particles.h5","r")
print("Keys: %s" % f.keys())
aeq2 = f["aeq2 2D"][()]
eta0 = f["initial_etas"][()]
time_sim = f["simulation time"][()]
lamda = f["lamda 2D"][()]
f.close();



norm = cm.colors.Normalize(vmin=eta0.min(), vmax=eta0.max())
cmap = cm.ScalarMappable(norm=norm, cmap=cm.jet)
cmap.set_array([])

fig, ax = plt.subplots()
for r in range(0,len(eta0)):
    ax.plot(lamda[r,:-91000]*R2D,aeq2[r,:-91000]*R2D,c=cmap.to_rgba(eta0[r]))

ax.grid(alpha=.3)
ax.set_xlabel('Latitude (deg)')
ax.set_ylabel('a_eq (deg)')
#ax.set_ylim(10,60)
#ax.set_xlim(19.5,26)
ticks=np.arange(0,2*np.pi,1)
cbar=fig.colorbar(cmap, ticks=ticks)
cbar.set_label('gyro-phase (rad)', rotation=270,labelpad=15)
ax.axvline(x=23,color="black", linestyle="--")
plt.show()




"""
dalpha=[]
for r in range(0,len(eta0)):
    ls=np.max(np.nonzero(aeq[r,:]))
    dalphaf=np.rad2deg(aeq[r,ls])-aeq0_deg
    dalpha.append(dalphaf)

fig, ax = plt.subplots()
ax.plot(eta0_deg,dalpha,color="tab:pink")
ax.grid(alpha=.3)
ax.set_xlim(0,360)
ax.set_xlabel('$ \eta_0$ (deg)')
ax.set_ylabel('$Da_{eq}$ (deg)')
ax.axvline(x=-5,color="black", linestyle="--")
plt.show()

fig, ax = plt.subplots()
for r in range(0,len(eta0)):
    ax.plot(lamda[r,:-91000]/RAD,E_kin[r,:-91000]/RAD,c=cmap.to_rgba(eta0[r]))
ax.grid(alpha=.3)
ax.set_xlim(-9,-0.1)
ax.set_xlabel('Latitude (deg)')
ax.set_ylabel('E [keV]')
ax.axvline(x=-5,color="black", linestyle="--")
plt.show()

fig, ax = plt.subplots()
for r in range(0,len(eta0)):
    ax.plot(lamda[r,:-91000]/RAD,mu_adiabatic[r,:-91000]/RAD,c=cmap.to_rgba(eta0[r]))
ax.grid(alpha=.3)
ax.set_xlim(-9,-0.1)
ax.set_ylim(6.585*10**(-47),6.60*10**(-47))
ax.set_xlabel('Latitude (deg)')
ax.set_ylabel('mu_adiabatic invariant')
ax.axvline(x=-5,color="black", linestyle="--")
plt.show()

fig, ax = plt.subplots()
for r in range(0,len(eta0)):
    ax.plot(lamda[r,:-91000]/RAD,upar[r,:-91000]/RAD,c=cmap.to_rgba(eta0[r]))
    ax.plot(lamda[r,:-91000]/RAD,vresz_o[r,:-91000]/RAD,color='tab:green')
ax.grid(alpha=.3)
ax.set_xlim(-9,-0.1)
ax.set_xlabel('Latitude (deg)')
ax.set_ylabel('v (m/s)')
ax.axvline(x=-5,color="black", linestyle="--")
plt.show()
"""