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
from scipy.stats import norm
import re

np.set_printoptions(threshold=sys.maxsize)

############################################# READ HDF5 ###################################################
#noWPI
f1 = h5py.File("h5files/detected_nowpi.h5","r")
#print("Keys: %s" % f1.keys())
detected_lamda = f1["ODPT.lamda"][()]
detected_time  = f1["ODPT.time"][()]
detected_id    = f1["ODPT.id"][()]
detected_alpha = f1["ODPT.alpha"][()]
telescope_lamda= f1["ODPT.latitude"][()]
population     = f1["population"][()]
t              = f1["t"][()]
#aeq0dstr       = f1["aeq0dstr"][()]
#da             = f1["da"][()]
lamda_start_d  = f1["lamda_start_d"][()]
lamda_end_d    = f1["lamda_end_d"][()]
aeq_start_d    = f1["aeq_start_d"][()]
aeq_end_d      = f1["aeq_end_d"][()]
Ekev0          = f1["Ekev0"][()]
By_wave        = f1["By_wave"][()]
aeq0           = f1["aeq0_plot"][()]
#alpha          = f1["alpha_plot"][()]
lamda0          = f1["lamda0_plot"][()]
#time           = f1["time_plot"][()]
#deta_dt        = f1["deta_dt"][()]
aeq0_bins       = f1["aeq0_bins"][()]

f1.close()

#noWPI and WPI afterwards
f2 = h5py.File("h5files/detected_both.h5","r")
#print("Keys: %s" % f2.keys())
detected_lamda_both = f2["ODPT.lamda"][()]
detected_time_both  = f2["ODPT.time"][()]
detected_id_both    = f2["ODPT.id"][()]
detected_alpha_both = f2["ODPT.alpha"][()]
telescope_lamda_both= f2["ODPT.latitude"][()]
population_both     = f2["population"][()]
t_both              = f2["t"][()]
#aeq0dstr_both       = f2["aeq0dstr"][()]
#da_both             = f2["da"][()]
lamda_start_d_both  = f2["lamda_start_d"][()]
lamda_end_d_both    = f2["lamda_end_d"][()]
aeq_start_d_both    = f2["aeq_start_d"][()]
aeq_end_d_both      = f2["aeq_end_d"][()]
Ekev0_both          = f2["Ekev0"][()]
By_wave_both        = f2["By_wave"][()]
aeq0_both           = f2["aeq0_plot"][()]
#alpha_both          = f2["alpha_plot"][()]
lamda0_both          = f2["lamda0_plot"][()]
#time_both           = f2["time_plot"][()]
#deta_dt_both        = f2["deta_dt"][()]
aeq0_bins_both       = f2["aeq0_bins"][()]

f2.close()


D2R=np.pi/180
R2D=1/D2R


##########################################################################################################
##########################################################################################################
###################################### POST PROCESSING - PLOTS ###########################################
##########################################################################################################
##########################################################################################################




############################# TELESCOPE SPECIFICATION && VARIABLES #######################################
time_bin  = 0.2                 #seconds to distinquish events(time resolution)
timesteps = int (t / time_bin)
view = 180 
sector_range = 15
sectors = int(view/sector_range)
########################################### FONTS AND COLORS #############################################
font = {'family': 'serif',
        'color':  'blue',
        'weight': 'bold',
        'size': 12}
colors = []
for i in range(max(timesteps,sectors)):      #colors to seperate timesteps or sectors.
    colors.append('#%06X' % randint(0, 0xFFFFFF))




######################################## PLOT INITIAL DISTRIBUTION #######################################
#aeq0 dstr
fig, ax = plt.subplots()
for sec in range(0,sectors):
    ax.scatter(sec,aeq0_bins[sec],s=2,alpha=1)
ax.grid(alpha=.3)
ax.set(xlabel="Sectors",xlim=(0,11),xticks=np.arange(0,12,1),ylabel="dN",title="Aeq0 distribution, sector range "+str(sector_range)+" degrees")
fig.savefig("simulation_MM/aeq0.png",dpi=200)

#aeq0-lamda0 after
fig, ax = plt.subplots()
ax.scatter(lamda0*R2D,aeq0*R2D,s=0.5,alpha=0.1)
ax.grid(alpha=.3)
ax.set(xlabel="Latitude(deg)",ylabel="Equatorial P.A",title="Initial lat-aeq of simulated particles",ylim=(1,179),xlim=(-90,90),xticks=np.linspace(-90,90,5))
ax.axhline(y = 90, color ="b", linestyle="dashed")
fig.savefig("simulation_MM/aeq0_lamda0.png",dpi=200)

################################### CROSSING PARTICLES LAMDA-TIME PLOT ####################################

#noWPI
fig, ax = plt.subplots()
ax.scatter(detected_time, detected_lamda*R2D, c = detected_id, s=0.3, cmap="viridis")
ax.grid(alpha=0.8)
ax.set(xlabel="time(s)", ylabel="latitude(deg)")
plt.title("$Population$: " +str(population)+ ", $lamda$: [" +str(lamda_start_d)+", "+str(lamda_end_d)+ "], $aeq0$: ["+str(aeq_start_d)+", "+str(aeq_end_d)+"]\nNorthward particles are captured below the satellite.\nSouthward particles are captured above the satellite",size="medium")
plt.annotate("SATELLITE",xy=(t/2,telescope_lamda+0.0002),color="blue",weight="semibold")
ax.ticklabel_format(useOffset=False)    #disable e notation.
ax.axhline(y = telescope_lamda ,color="b", linestyle="dashed")
plt.savefig("simulation_MM/Crossing_particles.png", dpi=100)

#BOTH
fig, ax = plt.subplots()
ax.scatter(detected_time_both, detected_lamda_both*R2D, c = detected_id_both, s=0.3, cmap="viridis")
ax.grid(alpha=0.8)
ax.set(xlabel="time(s)", ylabel="latitude(deg)")
plt.title("$Population$: " +str(population_both)+ ", $lamda$: [" +str(lamda_start_d_both)+", "+str(lamda_end_d_both)+ "], $aeq0$: ["+str(aeq_start_d_both)+", "+str(aeq_end_d_both)+"]\nNorthward particles are captured below the satellite.\nSouthward particles are captured above the satellite",size="medium")
plt.annotate("SATELLITE",xy=(t_both/2,telescope_lamda_both+0.0002),color="blue",weight="semibold")
ax.ticklabel_format(useOffset=False)    #disable e notation.
ax.axhline(y = telescope_lamda_both ,color="b", linestyle="dashed")
plt.savefig("simulation_MM/Crossing_particles_both.png", dpi=100)



############################################## BINNING ####################################################
sctr_flux = [ [0 for i in range(timesteps)] for j in range(sectors) ]   #sctr_flux[sectors][timesteps]
sum_flux = [0 for i in range(timesteps)]
for time,pa in zip(detected_time,detected_alpha): #Iterate in both array elements. No sorting required.
    timestep = math.floor(time/time_bin)
    sector   = math.floor(pa*R2D/sector_range)
    if(pa*R2D==180):
        sector = sectors-1 #to include p.a 180 in the last sector(11). Is this needed?
    sctr_flux[sector][timestep] += 1              #Number of detected particles in this sector-timestep.
    sum_flux[timestep] += 1

#BOTH
sctr_flux_both = [ [0 for i in range(timesteps)] for j in range(sectors) ]   #sctr_flux[sectors][timesteps]
sum_flux_both = [0 for i in range(timesteps)]
for time,pa in zip(detected_time_both,detected_alpha_both): #Iterate in both array elements. No sorting required.
    timestep = math.floor(time/time_bin)
    sector   = math.floor((pa*R2D)%180/sector_range)
    if(pa*R2D==180):
        sector = sectors-1 #to include p.a 180 in the last sector(11). Is this needed?
    sctr_flux_both[sector][timestep] += 1              #Number of detected particles in this sector-timestep.
    sum_flux_both[timestep] += 1
######################################### PARTICLE SUM - 360 PLOT ###########################################
fig, ax = plt.subplots()
plt.title("Detected particle sum in all look_dirs for "+str(t)+" seconds, in "+str(timesteps)+" timesteps\n Satellite @"+str(telescope_lamda)+" deg")
ax.set(xlabel="Time(s), in time_bins of "+str(time_bin)+"(s)", ylabel="Total Flux")
#ax.set_xticks(ticks=np.arange(0,t+time_bin,time_bin),minor=True) #ticks for time_bin seperation
ax.xaxis.grid(True, which='both')
for timestep in range(0,timesteps):
    if sum_flux[timestep]!=0:
        ax.scatter(timestep*time_bin+time_bin/2,sum_flux[timestep],s=10)
    if sum_flux_both[timestep]!=0:  
        ax.scatter(timestep*time_bin+time_bin/2,sum_flux_both[timestep],s=10, alpha=0.2)  #WPI would be transparent dots
plt.savefig("simulation_MM/Particle_sum.png", dpi=100)
###################################### (FLUX-P.A)*TIMESTEPS MOVIE ##########################################
fig,ax = plt.subplots()
FFMpegWriter = manimation.writers["ffmpeg"]
metadata2 = dict(title="P.A binning",comment="P.A bins of"+str(sector_range)+"degrees")
fps = 1
writer = FFMpegWriter(fps=fps, metadata = metadata2)
print("Generating P.A binning mp4 file...\nDuration of mp4 file will be:",(timesteps*fps), "seconds")
with writer.saving(fig, "simulation_MM/PA_binning.mp4", 100):
    for timestep in range(0,timesteps):           
            
        for sector in range(0,sectors): 
                        
            if sctr_flux[sector][timestep]!=0: #to show only non zero dots.
                ax.scatter(sector+0.5, sctr_flux[sector][timestep],c="black")  
                #If it differs with WPI, plot with red dots
                if (sctr_flux_both[sector][timestep] != sctr_flux[sector][timestep]): 
                    ax.scatter(sector+0.5, sctr_flux_both[sector][timestep],c="red")  

        ax.set_xticks(ticks=np.arange(0.5,sectors)) 
        ax.set_xticklabels(labels=np.arange(0,sectors),color="red",size="small")
        ax.set_xticks(ticks=np.arange(0,sectors),minor=True) 
        ax.set_xticklabels(labels=np.arange(0,view,sector_range),minor=True,size="small") 
        ax.set_ylim(0,np.amax(sctr_flux)+(0.1*np.amax(sctr_flux)) )
        ax.set_xlim(0,sectors)
        ax.set(ylabel="count")
        ax.set_title("P.A dstr, $time: "+str("{:.1f}".format(timestep*time_bin))+"s$", loc="left", size="small",color="blue",x=-0.15)
        ax.xaxis.grid(True, which='minor')

        writer.grab_frame()
        ax.clear() #clear data 








"""
######################################### LAST TIMESTEP - ALTER BINS ############################################
fig,ax = plt.subplots()

alter_time_bin  = [0.05,0.1,0.2]  #Try these time bins    
alter_sector_range = [7.5,15,30]  #With these sector ranges

FFMpegWriter = manimation.writers["ffmpeg"]
metadata2 = dict(title="Alter binning")
fps = 0.25
writer = FFMpegWriter(fps=fps, metadata = metadata2)
print("Generating alter binning mp4 file...\nDuration of mp4 file will be:",(len(alter_time_bin)*len(alter_sector_range)/fps), "seconds")

with writer.saving(fig, "simulation_MM/Alter_bins.mp4", 100):
   
    for srange in alter_sector_range:
        sectors = int(view/srange)
        for tbin in alter_time_bin:  
            timesteps = int (t / tbin)
            altered_sctr_flux = [ [0 for i in range(timesteps)] for j in range(sectors) ]   #altered_sctr_flux[sectors][timesteps]
            
            #Binning here
            for time,pa in zip(detected_time,detected_alpha): #Iterate in both array elements. No sorting required.
                timestep = math.floor(time/tbin)
                sector   = math.floor(pa*R2D/srange)
                if(pa*R2D==180):
                    sector = sectors-1 #to include p.a 180 in the last sector(11). Is this needed?
                altered_sctr_flux[sector][timestep] += 1              #Number of detected particles in this sector-timestep.
            #Plot here
            timestep = -1
            for sector in range(0,sectors): 
                if altered_sctr_flux[sector][timestep]!=0: #to show only non zero dots.
                    ax.scatter(sector+0.5, altered_sctr_flux[sector][timestep])  

            ax.set_xticks(ticks=np.arange(0.5,sectors)) 
            ax.set_xticklabels(labels=np.arange(0,sectors),color="red",size="small")
            ax.set_xticks(ticks=np.arange(0,sectors),minor=True) 
            if (srange>=15):            
                ax.set_xticklabels(labels=np.arange(0,view,srange),minor=True,size="small") 
                ax.set_xticklabels(labels=np.arange(0,view,srange),minor=True,size="small") 
            #ax.set_ylim(0,population/4)
            #ax.set_xlim(0,sectors)
            ax.set_title("Ending P.A dstr when\n$time$_$bin$: "+str(tbin)+"s, $sector$_$range$: "+str(srange)+"degrees", loc="left", size="small",color="blue",x=-0.15)
            ax.xaxis.grid(True, which='minor')
        
            writer.grab_frame()
            ax.clear() #clear data 
"""