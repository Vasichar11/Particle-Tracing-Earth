import numpy as np 
import matplotlib.pyplot as plt
from matplotlib import cm
import h5py


D2R=np.pi/180
R2D=1/D2R


#Reading the data from the file and then closing it.
f = h5py.File("./h5files/particles.h5","r")
print("Keys: %s" % f.keys())
aeqsu = f["aeqsu 2D"][()]
eta0 = f["initial_etas"][()]
time_sim = f["simulation time"][()]
lamda = f["lamda 2D"][()]
f.close();


######################### PLOT AEQSU-TIME MULTIPLE ETA ###############################
norm = cm.colors.Normalize(vmin=eta0.min(), vmax=eta0.max())
cmap = cm.ScalarMappable(norm=norm, cmap=cm.jet)
cmap.set_array([])
fig, ax = plt.subplots()
s=5
for r in range(0,len(eta0)):
    ax.plot(time_sim[:-2],aeqsu[r,:-2]*R2D,c=cmap.to_rgba(eta0[r]))
ax.grid(alpha=.3)
ax.set_xlabel('time (sec)')
ax.set_ylabel('a_eq (deg)')
ax.set_xlim(0.46,0.4775)
ax.set_ylim(69.95,70.15)
ticks=np.arange(0,2*np.pi,1)
cbar=fig.colorbar(cmap, ticks=ticks)
cbar.set_label('gyro-phase (rad)', rotation=270,labelpad=15)
fig.savefig('./multi_plots_cc/aeq_time_h5.png')
plt.show()
######################### PLOT AEQSU-TIME MULTIPLE ETA ################################
"""
######################### PLOT AEQSU-LAMDA MULTIPLE ETA ################################
norm = cm.colors.Normalize(vmin=eta0.min(), vmax=eta0.max())
cmap = cm.ScalarMappable(norm=norm, cmap=cm.jet)
cmap.set_array([])
fig, ax = plt.subplots()
for r in range(0,len(eta0)):
    ax.plot(lamda[r,:-2]*R2D,aeqsu[r,:-2]*R2D,c=cmap.to_rgba(eta0[r]))
ax.grid(alpha=.3)
#ax.set_xlim(-9,-2)
ax.set_xlabel('Latitude (deg)')
ax.set_ylabel('a_eq (deg)')
#ax.set_ylim(69.96,70.04)
cbar=fig.colorbar(cmap, ticks=ticks)
cbar.set_label('gyro-phase (rad)', rotation=270,labelpad=15)
ax.axvline(x=-5,color="black", linestyle="--")
fig.savefig('./multi_plots_cc/aeq_lamda_h5.png')
plt.show()
######################### PLOT AEQSU-LAMDA MULTIPLE ETA: DONE ################################
"""
