import numpy as np 
import matplotlib.pyplot as plt
from matplotlib import cm
import h5py

D2R=np.pi/180
R2D=1/D2R

#Reading the data from the file and then closing it.
f = h5py.File("particles.h5","r")
print("Keys: %s" % f.keys())
aeq2 = f["aeq2 2D"][()]
aeq = f["aeq 2D"][()]
eta0 = f["initial_etas"][()]
time_sim = f["simulation time"][()]
lamda = f["lamda 2D"][()]
aeq0_deg = f["aeq0 single scalar"][()]
f.close();
eta0_deg = eta0*R2D

################################# LAMDA-AEQ-ETA PLOT ##################################
norm = cm.colors.Normalize(vmin=eta0.min(), vmax=eta0.max())
cmap = cm.ScalarMappable(norm=norm, cmap=cm.jet)
cmap.set_array([])
fig, ax = plt.subplots()
for r in range(0,len(eta0)):
    ax.plot(lamda[r,:-91000]*R2D,aeq2[r,:-91000]*R2D,c=cmap.to_rgba(eta0[r]))
ax.grid(alpha=.3)
ax.set_xlabel('Latitude (deg)')
ax.set_ylabel('a_eq (deg)')

ax.set_xlim(-30,-0.1)
ax.set_ylim(10,60)

ticks=np.arange(2.88,3.4,0.04)
#ticks=np.arange(eta0[0],eta0[-1],eta0[1]-eta0[0]) #arrange from first to last eta
cbar=fig.colorbar(cmap, ticks=ticks)
cbar.set_label('gyro-phase (rad)', rotation=270,labelpad=15)
ax.axvline(x=-5,color="black", linestyle="--")
plt.show()
################################# LAMDA-AEQ-ETA PLOT ##################################

###################################### ΔA-ETA PLOT#####################################
dalpha=[]
for r in range(0,len(eta0)):
    ls=np.max(np.nonzero(aeq[r,:]))#Return the indices of the elements that are non-zero and Return max
    dalphaf=np.rad2deg(aeq[r,ls])-aeq0_deg
    dalpha.append(dalphaf)#Δa
fig, ax = plt.subplots()
ax.plot(eta0_deg,dalpha,color="tab:pink")
# ax.legend()
ax.grid(alpha=.3)

ax.set_xlim(175,185) 
# ax.set_ylim(168.27,168.32)

ax.set_xlabel('$ \eta_0$ (deg)')
ax.set_ylabel('$Da_{eq}$ (deg)')
ax.axvline(x=-5,color="black", linestyle="--")
plt.show()
##################################### ΔA-ETA PLOT######################################


#Max dalpha
print(eta0_deg[np.argmax(dalpha)])
