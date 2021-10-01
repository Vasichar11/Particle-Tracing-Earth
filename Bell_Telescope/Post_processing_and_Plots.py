import numpy as np 
import random
import sys
import matplotlib.pyplot as plt
from matplotlib import cm
import h5py
import math
import matplotlib.animation as manimation
from random import randint
from matplotlib.ticker import MaxNLocator
from collections import Counter
import math as mt

D2R=np.pi/180
R2D=1/D2R

np.set_printoptions(threshold=sys.maxsize)


############################################# READ DATA ###################################################

#noWPI
f1 = h5py.File("h5files/detected.h5","r")
#print("Keys: %s" % f1.keys())
detected_lamda = f1["ODPT.lamda"][()]
detected_time  = f1["ODPT.time"][()]
detected_id    = f1["ODPT.id"][()]
detected_alpha = f1["ODPT.alpha"][()]
telescope_lamda= f1["ODPT.latitude"][()]
population     = f1["population"][()]
t              = f1["t"][()]
lamda_start_d  = f1["lamda_start_d"][()]
lamda_end_d    = f1["lamda_end_d"][()]
aeq_start_d    = f1["aeq_start_d"][()]
aeq_end_d      = f1["aeq_end_d"][()]
Ekev0          = f1["Ekev0"][()]
By_wave        = f1["By_wave"][()]
#aeq            = f1["aeq_plot"][()]
#alpha          = f1["alpha_plot"][()]
#lamda          = f1["lamda_plot"][()]
#time           = f1["time_plot"][()]
#deta_dt        = f1["deta_dt"][()]

f1.close();

"""
#WPI
f2 = h5py.File("h5files/detected_WPI.h5","r")
#print("Keys: %s" % f2.keys())
detected_lamda_WPI = f2["ODPT.lamda"][()]
detected_time_WPI  = f2["ODPT.time"][()]
detected_id_WPI    = f2["ODPT.id"][()]
detected_alpha     = f2["ODPT.alpha"][()]
telescope_lamda_WPI= f2["ODPT.latitude"][()]
population_WPI     = f2["population"][()]
t_WPI              = f2["t"][()]
lamda_start_d_WPI  = f1["lamda_start_d"][()]
lamda_end_d_WPI    = f1["lamda_end_d"][()]
aeq_start_d_WPI    = f1["aeq_start_d"][()]
aeq_end_d_WPI      = f1["aeq_end_d"][()]
Ekev0_WPI          = f1["Ekev0"][()]
By_wave_WPI        = f1["By_wave"][()]
aeq_WPI            = f2["aeq_plot"][()]
alpha_WPI          = f2["alpha_plot"][()]
lamda_WPI          = f2["lamda_plot"][()]
time_WPI           = f2["time_plot"][()]
deta_dt_WPI        = f2["deta_dt"][()]

f2.close();
"""

##########################################################################################################
###################################### POST PROCESSING - PLOTS ###########################################
##########################################################################################################

########################################## COMPARE WPI - noWPI ###########################################
"""
detected_id     = list(detected_id)     #Turn into lists
detected_id_WPI = list(detected_id_WPI)
new_particles = list(set(detected_id_WPI) - set(detected_id))  #Particles that were detected with WPI but not before.
print("New id's when WPI",By_wave_WPI,"\n",new_particles)
print(new_particles)
new_part = new_particles[1] #take one of them

#most = max(set(detected_id_WPI), key = detected_id_WPI.count)  #One of the most detected particles(regularly bouncing whole time).

fig, ax = plt.subplots(3,2)  

#noWPI
ax[0,0].set_title("noWPI",size="small",color="darkorange")     #AEQ-TIME
ax[0,0].plot(time[new_part,:-1],aeq[new_part,:-1]*R2D)
ax[0,0].grid(alpha=.3)
ax[0,0].set_ylim(10,80)
ax[0,0].set_ylabel('$aeq (deg)$',size="small")

ax[1,0].plot(time[new_part,1:-1],deta_dt[new_part,1:-1])        #DETA_DT-TIME
ax[1,0].grid(alpha=.3)
ax[1,0].set_ylabel('$deta_dt$',size="small")

ax[2,0].plot(time[new_part,:-1],lamda[new_part,:-1]*R2D)             #LAMDA-TIME
ax[2,0].axhline(y = telescope_lamda, color ="b", linestyle="dashed")
ax[2,0].set_ylim(-40,40)
ax[2,0].annotate("SATELLITE",xy=(0,telescope_lamda_WPI+0.2),color="blue",weight="semibold",size="small")
ax[2,0].grid(alpha=.3)
ax[2,0].set_ylabel('$lamda (deg)$',size="small")
ax[2,0].set_xlabel('$time (sec)$',size="x-small")


#WPI
ax[0,1].set_title("WPI",size="small",color="darkorange")         #AEQ-TIME
ax[0,1].plot(time_WPI[new_part,:-1],aeq_WPI[new_part,:-1]*R2D)
ax[0,1].set_ylim(10,80)
ax[0,1].grid(alpha=.3)

ax[1,1].plot(time_WPI[new_part,1:-1],deta_dt_WPI[new_part,1:-1])    #DETA_DT-TIME
ax[1,1].grid(alpha=.3)

ax[2,1].plot(time_WPI[new_part,:-1],lamda_WPI[new_part,:-1]*R2D)
ax[2,1].axhline(y = telescope_lamda, color ="b", linestyle="dashed")     #LAMDA-TIME
ax[2,1].set_ylim(-40,40)
ax[2,1].annotate("SATELLITE",xy=(0,telescope_lamda_WPI+0.8),color="blue",weight="semibold",size="small")
ax[2,1].grid(alpha=.3)
ax[2,1].set_xlabel('$time (sec)$',size="x-small")

fig.savefig('simulation_MM/aeq_deta_lamda.png' , dpi=200, facecolor='white',transparent=False)

"""
################################### CROSSING PARTICLES LAMDA-TIME PLOT ####################################

#noWPI
fig, ax = plt.subplots()
ax.scatter(detected_time, detected_lamda*R2D, c = detected_id, s=0.3, cmap="viridis")
ax.grid(alpha=0.8)
ax.set(xlabel="time(s)", ylabel="latitude(deg)")
plt.title("$Population$: " +str(population)+ ", $lamda$: [" +str(lamda_start_d)+", "+str(lamda_end_d)+ "], $aeq0$: ["+str(aeq_start_d)+", "+str(aeq_end_d)+"], $By wave$: "+str(By_wave)+"\nNorthward particles are captured below the satellite.\nSouthward particles are captured above the satellite",size="medium")
plt.annotate("SATELLITE",xy=(t/2,telescope_lamda+0.0002),color="blue",weight="semibold")
ax.ticklabel_format(useOffset=False)    #disable e notation.
ax.axhline(y = telescope_lamda ,color="b", linestyle="dashed")
plt.savefig("simulation_MM/Crossing_particles.png", dpi=100)
"""

#WPI
fig, ax = plt.subplots()
ax.scatter(detected_time_WPI, detected_lamda_WPI*R2D, c = detected_id_WPI, s=7, cmap="viridis")
ax.grid(alpha=0.8)
ax.set(xlabel="time(s)", ylabel="latitude(deg)")
plt.title("Population: " +str(population_WPI)+ ", lamda: [" +str(lamda_start_d_WPI)+", "+str(lamda_end_d_WPI)+ "], aeq0: ["+str(aeq_start_d_WPI)+", "+str(aeq_end_d_WPI)+"], "+str(By_wave_WPI)+"\nNorthward particles are captured below the satellite.\nSouthward particles are captured above the satellite",size="medium")
plt.annotate("SATELLITE",xy=(t/2,telescope_lamda_WPI+0.0002),color="blue",weight="semibold")
ax.ticklabel_format(useOffset=False)    #disable e notation.
ax.axhline(y = telescope_lamda_WPI ,color="b", linestyle="dashed")
plt.savefig("simulation_MM/Crossing_particles_WPI.png", dpi=100)

#PLOT ONLY NEW PARTICLES #NOT VALID...
fig, ax = plt.subplots()
j=0
for id in new_particles:
    ax.scatter(detected_time_WPI[id], detected_lamda_WPI[id]*R2D, c = "r", s=7)
    j+=1
ax.grid(alpha=0.8)
ax.set(xlabel="time(s)", ylabel="latitude(deg)")
plt.title("New particles : " +str(len(new_particles))+ "\n"+"Northward particles are captured below the satellite.\nSouthward particles are captured above the satellite.",size="medium")
plt.annotate("SATELLITE",xy=(t/2,telescope_lamda_WPI+0.0002),color="blue",weight="semibold")
ax.ticklabel_format(useOffset=False)    #disable e notation.
ax.axhline(y = telescope_lamda_WPI ,color="b", linestyle="dashed")

plt.savefig("simulation_MM/Crossing_particles_WPI_onlyNEW.png", dpi=100)

"""
####################################### TELESCOPE SPECIFICATION ##########################################

time_bin  = 0.1                 #seconds to distinquish events(time resolution)
timesteps = int (t / time_bin)
view = 180 
sector_range = 15
sectors = int(view/sector_range)

############################################## BINNING ####################################################


sctr_flux = [ [0 for i in range(timesteps)] for j in range(sectors) ]   #sctr_flux[sectors][timesteps]
sum_flux = [0 for i in range(timesteps)]

for time,pa in zip(detected_time,detected_alpha): #Iterate in both array elements. No sorting required.
    timestep = math.floor(time*10)
    sector   = math.floor(pa*R2D/sector_range)
    if(pa==180):
        sector = sectors-1 #to include p.a 180 in the last sector(11). Is this needed?
    sctr_flux[sector][timestep] += 1              #Number of detected particles in this sector-time_bin.
    sum_flux[timestep] += 1



######################################### PARTICLE SUM - 360 PLOT ###########################################

fig, ax = plt.subplots()
plt.title("Detected particle sum in all look_dirs for "+str(t)+" seconds, in "+str(timesteps)+" timesteps\n Satellite @"+str(telescope_lamda)+" deg")
ax.set(xlabel="Time(s), in time_bins of "+str(time_bin)+"(s)", ylabel="Total Flux")
ax.set_xticks(ticks=np.arange(0,t+time_bin,time_bin),minor=True) #ticks for time_bin seperation
ax.xaxis.grid(True, which='both')
for timestep in range(0,timesteps):
    if sum_flux[timestep]!=0:
        ax.scatter(timestep*0.1+0.1/2,sum_flux[timestep],s=10)

plt.savefig("simulation_MM/Particle_sum.png", dpi=100)

####################################### (FLUX-TIME)*SECTORS  MOVIE ##########################################

font = {'family': 'serif',
        'color':  'blue',
        'weight': 'medium',
        'size': 10,
        }
colors = []
for i in range(max(timesteps,sectors)):      #colors to seperate timesteps or sectors.
    colors.append('#%06X' % randint(0, 0xFFFFFF))

fig,ax = plt.subplots()
FFMpegWriter = manimation.writers["ffmpeg"]
metadata = dict(title="Time binning",comment="Time bins of 0.1s")
fps = 1
writer = FFMpegWriter(fps=fps, metadata = metadata)
print("Generating Time binning mp4 file...\nDuration of mp4 file will be:",(sectors*fps), "seconds")


with writer.saving(fig, "simulation_MM/Time_binning.mp4", 100):

    for sector in range(0,sectors):

        for timestep in range(0,timesteps):           
                
            if sctr_flux[sector][timestep] !=0: #to show only non zero dots.
                ax.scatter(timestep*time_bin+(time_bin/2),sctr_flux[sector][timestep],c=colors[timestep],s=10)
       
        ax.set(xlabel="Time(s), in time_bins of "+str(time_bin)+" seconds", ylabel="Sector Flux")
        ax.set_xticks(ticks=np.arange(0,t+time_bin,time_bin),minor=True) #ticks for time_bin seperation
        ax.xaxis.grid(True, which='both')
        ax.set_ylim(0,np.amax(sctr_flux))
        ax.set_title("Particle Flux in Sector "+str(sector)+" - Look direction: ("+str(180+sector*sector_range)+"\N{DEGREE SIGN},"+str(180+sector*sector_range+sector_range)+"\N{DEGREE SIGN}) Simulation time: "+str(t)+"s\n                                         -Detectable P.A: [ "+str(sector*sector_range)+"\N{DEGREE SIGN},"+str(sector*sector_range+sector_range)+"\N{DEGREE SIGN}) Satellite @"+str(telescope_lamda)+" deg",loc="left",fontdict=font)              


        writer.grab_frame()
        ax.clear() #clear data 










