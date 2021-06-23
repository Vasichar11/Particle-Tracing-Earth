import numpy as np 
import matplotlib.pyplot as plt
from matplotlib import cm

eta0 = np.loadtxt("./cc_text/eta0.txt")
time_sim = np.loadtxt("./cc_text/time_sim.txt")

part = np.loadtxt("./cc_text/particle.txt")
#Columns(see outfile): lamda->0-population, aeq->population-2*population, aeqsu->2*population-3*population, alpha->3*population-4*population 


varsOfpart0 = np.loadtxt("./cc_text/varsOfpart0.txt")
varsOfpart1 = np.loadtxt("./cc_text/varsOfpart1.txt")
varsOfpart2 = np.loadtxt("./cc_text/varsOfpart2.txt")
varsOfpart3 = np.loadtxt("./cc_text/varsOfpart3.txt")
varsOfpart4 = np.loadtxt("./cc_text/varsOfpart4.txt")
varsOfpart5 = np.loadtxt("./cc_text/varsOfpart5.txt")
#Columns(see runge kutta): mu_ad_li->0 , Ftheta_o->1, Bxw_out->2, Byw_out->3, Bzw_out->4, gamma_out->5.

D2R=np.pi/180
R2D=1/D2R

################################## MU PLOT #####################################################
#In[40]:
norm = cm.colors.Normalize(vmin=eta0.min(), vmax=eta0.max())
cmap = cm.ScalarMappable(norm=norm, cmap=cm.jet)
cmap.set_array([])
fig, ax = plt.subplots()
s=5
ax.plot(time_sim[:-2],varsOfpart0[:-2,0],c=cmap.to_rgba(eta0[0]))
ax.plot(time_sim[:-2],varsOfpart1[:-2,0],c=cmap.to_rgba(eta0[1]))
ax.plot(time_sim[:-2],varsOfpart2[:-2,0],c=cmap.to_rgba(eta0[2]))
ax.plot(time_sim[:-2],varsOfpart3[:-2,0],c=cmap.to_rgba(eta0[3]))
ax.plot(time_sim[:-2],varsOfpart4[:-2,0],c=cmap.to_rgba(eta0[4]))
ax.plot(time_sim[:-2],varsOfpart5[:-2,0],c=cmap.to_rgba(eta0[5]))
ax.grid(alpha=.3)
ax.set_xlim(0.46,0.49)
ax.set_xlabel('time (sec)')
ax.set_ylabel('$\mu$ ')
ticks=np.arange(0,2*np.pi,1)
cbar=fig.colorbar(cmap, ticks=ticks)
fig.savefig('./multi_plots_cc/mu.png')
################################## MU PLOT #####################################################



################################## FTHETA PLOT #####################################################
fig, ax = plt.subplots()
s=5
ax.plot(time_sim[:-2],varsOfpart0[:-2,1],c=cmap.to_rgba(eta0[0]))
ax.plot(time_sim[:-2],varsOfpart1[:-2,1],c=cmap.to_rgba(eta0[1]))
ax.plot(time_sim[:-2],varsOfpart2[:-2,1],c=cmap.to_rgba(eta0[2]))
ax.plot(time_sim[:-2],varsOfpart3[:-2,1],c=cmap.to_rgba(eta0[3]))
ax.plot(time_sim[:-2],varsOfpart4[:-2,1],c=cmap.to_rgba(eta0[4]))
ax.plot(time_sim[:-2],varsOfpart5[:-2,1],c=cmap.to_rgba(eta0[5]))
ax.grid(alpha=.3)
ax.set_xlabel('time (sec)')
ax.set_ylabel('Ftheta ')
ticks=np.arange(0,2*np.pi,1)
cbar=fig.colorbar(cmap, ticks=ticks)
cbar.set_label('gyro-phase (rad)', rotation=270,labelpad=15)
fig.savefig('./multi_plots_cc/Ftheta.png')
################################## FTHETA PLOT #####################################################



################################## AEQ TIME PLOT #####################################################
population = 6;
norm = cm.colors.Normalize(vmin=eta0.min(), vmax=eta0.max())
cmap = cm.ScalarMappable(norm=norm, cmap=cm.jet)
cmap.set_array([])
fig, ax = plt.subplots()
s=5
ax.plot(time_sim[:-2],part[:-2, population + 0]*R2D,c=cmap.to_rgba(eta0[0]))
ax.plot(time_sim[:-2],part[:-2, population + 1]*R2D,c=cmap.to_rgba(eta0[1]))
ax.plot(time_sim[:-2],part[:-2, population + 2]*R2D,c=cmap.to_rgba(eta0[2]))
ax.plot(time_sim[:-2],part[:-2, population + 3]*R2D,c=cmap.to_rgba(eta0[3]))
ax.plot(time_sim[:-2],part[:-2, population + 4]*R2D,c=cmap.to_rgba(eta0[4]))
ax.plot(time_sim[:-2],part[:-2, population + 5]*R2D,c=cmap.to_rgba(eta0[5]))
ax.grid(alpha=.3)
ax.set_xlim(0,3.1)
ax.set_xlabel('time (sec)')
ax.set_ylabel('a_eq (deg)')
ax.set_ylim(69.98,70.02)
ticks=np.arange(0,2*np.pi,1)
cbar=fig.colorbar(cmap, ticks=ticks)
cbar.set_label('gyro-phase (rad)', rotation=270,labelpad=15)
fig.savefig('./multi_plots_cc/aeq_time.png')
################################## AEQ TIME PLOT #####################################################



################################## ENERGY PLOT #####################################################
m_e=9.10938291e-31   #electron mass in kg
c=2.997956376932163e8 #speed of light in m/s

norm = cm.colors.Normalize(vmin=eta0.min(), vmax=eta0.max())
cmap = cm.ScalarMappable(norm=norm, cmap=cm.jet)
cmap.set_array([])
fig, ax = plt.subplots()
s=5				#Energy=m_e*c*c* (gamma_out          -1 ) *6.2415*(10**15) 
ax.plot(time_sim[:-2],  m_e*c*c* (varsOfpart0[:-2,5] -1 ) *6.2415*(10**15)   ,c=cmap.to_rgba(eta0[0]))
ax.plot(time_sim[:-2],  m_e*c*c* (varsOfpart1[:-2,5] -1 ) *6.2415*(10**15)   ,c=cmap.to_rgba(eta0[1]))
ax.plot(time_sim[:-2],  m_e*c*c* (varsOfpart2[:-2,5] -1 ) *6.2415*(10**15)   ,c=cmap.to_rgba(eta0[2]))
ax.plot(time_sim[:-2],  m_e*c*c* (varsOfpart3[:-2,5] -1 ) *6.2415*(10**15)   ,c=cmap.to_rgba(eta0[3]))
ax.plot(time_sim[:-2],  m_e*c*c* (varsOfpart4[:-2,5] -1 ) *6.2415*(10**15)   ,c=cmap.to_rgba(eta0[4]))
ax.plot(time_sim[:-2],  m_e*c*c* (varsOfpart5[:-2,5] -1 ) *6.2415*(10**15)   ,c=cmap.to_rgba(eta0[5]))
ax.grid(alpha=.3)
ax.set_xlim(0,1)
ax.set_xlabel('time (sec)')
ax.set_ylabel('Energy (keV)')
ticks=np.arange(0,2*np.pi,1)
cbar=fig.colorbar(cmap, ticks=ticks)
cbar.set_label('gyro-phase (rad)', rotation=270,labelpad=15)
fig.savefig('./multi_plots_cc/Energy.png')
################################## ENERGY PLOT #####################################################


################################## SINGLE AEQ_SU PLOT #####################################################
fig, ax = plt.subplots()
s=5
ax.plot(time_sim[:-2],part[:-2,2*population+0]*R2D,c=cmap.to_rgba(eta0[0]))
ax.grid(alpha=.3)
ax.set_xlim(0,3.1)
ax.set_xlabel('time (sec)')
ax.set_ylabel('$ a_{eq}$ (deg)')
ax.set_xlim(0.46,0.4775)
ax.set_ylim(69.95,70.15)
ticks=np.arange(0,2*np.pi,1)
cbar=fig.colorbar(cmap, ticks=ticks)
cbar.set_label('gyro-phase (rad)', rotation=270,labelpad=15)
fig.savefig('./multi_plots_cc/single_aeqsu_gyro.png')
################################## SINGLE AEQ_SU PLOT #####################################################



################################## ALPHA LAT PLOT #####################################################
# In[32]:
norm = cm.colors.Normalize(vmin=eta0.min(), vmax=eta0.max())
cmap = cm.ScalarMappable(norm=norm, cmap=cm.jet)
cmap.set_array([])
fig, ax = plt.subplots()
s=5
ax.plot(time_sim[:-2], part[:-2, 3*population +0]*R2D, c=cmap.to_rgba(eta0[0]))
ax.plot(time_sim[:-2], part[:-2, 3*population +1]*R2D, c=cmap.to_rgba(eta0[1]))
ax.plot(time_sim[:-2], part[:-2, 3*population +2]*R2D, c=cmap.to_rgba(eta0[2]))
ax.plot(time_sim[:-2], part[:-2, 3*population +3]*R2D, c=cmap.to_rgba(eta0[3]))
ax.plot(time_sim[:-2], part[:-2, 3*population +4]*R2D, c=cmap.to_rgba(eta0[4]))
ax.plot(time_sim[:-2], part[:-2, 3*population +5]*R2D, c=cmap.to_rgba(eta0[5]))
ax.grid(alpha=.3)
ax.set_xlabel('Latitude (deg)')
ax.set_ylabel('a_eq (deg)')
ax.set_xlim(0,5)
ax.set_ylim(109.5,110.4)
ticks=np.arange(0,2*np.pi,1)
cbar=fig.colorbar(cmap, ticks=ticks)
cbar.set_label('gyro-phase (rad)', rotation=270,labelpad=15)
fig.savefig('./multi_plots_cc/alpha_lat.png')
################################## ALPHA LAT PLOT #####################################################


################################## SECOND PARTICLE'S FIELDS PLOT #####################################################
# In[36]:
fig, ax = plt.subplots()
ax.plot(time_sim[:-2], varsOfpart1[:-2,2], label='Bxw')
ax.plot(time_sim[:-2], varsOfpart1[:-2,3], label='Byw')
ax.plot(time_sim[:-2], varsOfpart1[:-2,4], label='Bzw')
ax.legend()
ax.grid(alpha=.3)
ax.set_xlabel('$ time $ (s)')
ax.set_ylabel('Wave magnetic field (T) ')
fig.savefig('./multi_plots_cc/fields.png')
################################## SECOND PARTICLE'S FIELDS PLOT #####################################################





