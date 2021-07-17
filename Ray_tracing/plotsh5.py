import numpy as np 
import matplotlib.pyplot as plt
from matplotlib import cm
import h5py


D2R=np.pi/180
R2D=1/D2R
m_e=9.10938291e-31   #electron mass in kg
c=2.997956376932163e8 #speed of light in m/s



#Reading the data from the file and then closing it.
f = h5py.File("particles.h5","r")
print("Keys: %s" % f.keys())
aeqsu = f["aeqsu 2D"][()]
eta0 = f["initial_etas"][()]
time_sim = f["simulation time"][()]
lamda = f["lamda 2D"][()]
mu_ad_li = f["mu_ad_li 2D"][()]
alpha = f["alpha 2D"][()]
Ftheta_o = f["Ftheta 2D"][()]
gamma_out = f["gamma 2D"][()]
Bxw_out = f["Bx 2D"][()]
Byw_out = f["By 2D"][()]
Bzw_out = f["Bz 2D"][()]

f.close();


from matplotlib import cm
norm = cm.colors.Normalize(vmin=eta0.min(), vmax=eta0.max())
cmap = cm.ScalarMappable(norm=norm, cmap=cm.jet)
cmap.set_array([])

fig, ax = plt.subplots()
s=5

for r in range(0,len(eta0)):
    ax.plot(time_sim[r,:-2],mu_ad_li[r,:-2],c=cmap.to_rgba(eta0[r]))
ax.grid(alpha=.3)
ax.set_xlim(0.46,0.49)
ax.set_xlabel('time (sec)')
ax.set_ylabel('$\mu$ ')
ticks=np.arange(0,2*np.pi,1)
cbar=fig.colorbar(cmap, ticks=ticks)
cbar.set_label('gyro-phase (rad)', rotation=270,labelpad=15)
fig.savefig('mu.png' , dpi=200, facecolor='white',transparent=False)
plt.show()




fig, ax = plt.subplots()
s=5

for r in range(0,len(eta0)):
    ax.plot(time_sim[r,:-2],Ftheta_o[r,:-2],c=cmap.to_rgba(eta0[r]))
ax.grid(alpha=.3)
ax.set_xlabel('time (sec)')
ax.set_ylabel('Ftheta')
ticks=np.arange(0,2*np.pi,1)
cbar=fig.colorbar(cmap, ticks=ticks)
cbar.set_label('gyro-phase (rad)', rotation=270,labelpad=15)
plt.show()



""" Input h5 aeq.
norm = cm.colors.Normalize(vmin=eta0.min(), vmax=eta0.max())
cmap = cm.ScalarMappable(norm=norm, cmap=cm.jet)
cmap.set_array([])

fig, ax = plt.subplots()
s=5

for r in range(0,len(eta0)):
    ax.plot(time_sim[r,:-2],aeq[r,:-2]/D2R,c=cmap.to_rgba(eta0[r]))
ax.grid(alpha=.3)
ax.set_xlim(0,3.1)
ax.set_xlabel('time (sec)')
ax.set_ylabel('a_eq (deg)')
ax.set_ylim(69.98,70.02)
ticks=np.arange(0,2*np.pi,1)
cbar=fig.colorbar(cmap, ticks=ticks)
cbar.set_label('gyro-phase (rad)', rotation=270,labelpad=15)
plt.show()
"""


Energy=m_e*c*c*(gamma_out-1)*6.2415*(10**15)


norm = cm.colors.Normalize(vmin=eta0.min(), vmax=eta0.max())
cmap = cm.ScalarMappable(norm=norm, cmap=cm.jet)
cmap.set_array([])

fig, ax = plt.subplots()
s=5


for r in range(0,len(eta0)):
    ax.plot(time_sim[r,:-2],Energy[r,:-2],c=cmap.to_rgba(eta0[r]))
ax.grid(alpha=.3)
ax.set_xlim(0,1)
ax.set_xlabel('time (sec)')
ax.set_ylabel('Energy (keV)')
ticks=np.arange(0,2*np.pi,1)
cbar=fig.colorbar(cmap, ticks=ticks)
cbar.set_label('gyro-phase (rad)', rotation=270,labelpad=15)
fig.savefig('energy.png' , dpi=200, facecolor='white',transparent=False)
plt.show()

fig, ax = plt.subplots()
s=5


for r in range(0,len(eta0)-5):
    ax.plot(time_sim[r,:-2],aeqsu[r,:-2]/D2R,c=cmap.to_rgba(eta0[r]))
ax.grid(alpha=.3)
ax.set_xlim(0,3.1)
ax.set_xlabel('time (sec)')
ax.set_ylabel('$ a_{eq}$ (deg)')
ax.set_xlim(0.46,0.4775)
ax.set_ylim(69.95,70.15)
ticks=np.arange(0,2*np.pi,1)
cbar=fig.colorbar(cmap, ticks=ticks)
cbar.set_label('gyro-phase (rad)', rotation=270,labelpad=15)
fig.savefig('a_eq.png' , dpi=200, facecolor='white',transparent=False)
plt.show()


norm = cm.colors.Normalize(vmin=eta0.min(), vmax=eta0.max())
cmap = cm.ScalarMappable(norm=norm, cmap=cm.jet)
cmap.set_array([])

fig, ax = plt.subplots()
s=5


for r in range(0,len(eta0)):
    ax.plot(time_sim[r,:-2],alpha[r,:-2]/D2R,c=cmap.to_rgba(eta0[r]))
ax.grid(alpha=.3)
ax.set_xlabel('Latitude (deg)')
ax.set_ylabel('a_eq (deg)')
ax.set_xlim(0,5)
ax.set_ylim(109.5,110.4)
ticks=np.arange(0,2*np.pi,1)
cbar=fig.colorbar(cmap, ticks=ticks)
cbar.set_label('gyro-phase (rad)', rotation=270,labelpad=15)
plt.show()

fig, ax = plt.subplots()
ax.plot(time_sim[1,:-2],Bxw_out[1,:-2],label='Bxw')
ax.plot(time_sim[1,:-2],Byw_out[1,:-2],label='Byw')
ax.plot(time_sim[1,:-2],Bzw_out[1,:-2],label='Bzw')
ax.legend()
ax.grid(alpha=.3)
ax.set_xlabel('$ time $ (s)')
ax.set_ylabel('Wave magnetic field (T) ')
plt.show()





