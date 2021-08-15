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
eta = f["eta 2D"][()]
time_sim = f["simulation time"][()]
lamda = f["lamda 2D"][()]
population = f["population"][()]
f.close();

lamda0 = lamda[:,0]
eta0   = eta[:,0]
aeq0   = aeq[:,0]
eta0_deg = eta0*R2D
aeq0_deg = aeq0*R2D
lamda0_deg = lamda0*R2D


################################# LAMDA-AEQ-ETA PLOT ##################################
norm = cm.colors.Normalize(vmin=eta0_deg.min(), vmax=eta0_deg.max())
cmap = cm.ScalarMappable(norm=norm, cmap=cm.jet)
cmap.set_array([])
fig, ax = plt.subplots()
for r in range(0,population):
    ax.plot(lamda[r,:-91000]*R2D,aeq2[r,:-91000]*R2D,c=cmap.to_rgba(eta0_deg[r]))


#ticks=np.arange(2.88,3.4,0.04)
#ticks=np.arange(eta0[0],eta0[-1],eta0[1]-eta0[0]) #arrange from first to last eta
ticks=np.arange(eta0_deg.min(), eta0_deg.max() + 90 ,90) #add step(180) to end of interval to make it full-closed interval.
cbar=fig.colorbar(cmap, ticks=ticks)
cbar.set_label('eta0(deg)', rotation=270,labelpad=15)
ax.axvline(x=-25,color="black", linestyle="--")
ax.grid(alpha=.3)
ax.set_xlabel('Latitude (deg)')
ax.set_ylabel('a_eq (deg)')

#ax.set_xlim(-30,-20)
#ax.set_ylim(29.8,30.2)

plt.show()
################################# LAMDA-AEQ-ETA PLOT ##################################


