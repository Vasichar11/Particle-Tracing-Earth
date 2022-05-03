from ipaddress import summarize_address_range
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
import csv

np.set_printoptions(threshold=sys.maxsize)
D2R=np.pi/180
R2D=1/D2R
print(len(np.arange(0.5,90+0.5)))
############################################# READ HDF5 ###################################################
#noWPI read
f1 = h5py.File("h5files/20000p_no_wpi.h5","r")
#print("Keys: %s" % f1.keys())
detected_lamda = f1["ODPT.lamda"][()]
detected_time  = f1["ODPT.time"][()]
detected_id    = f1["ODPT.id"][()]
detected_alpha = f1["ODPT.alpha"][()]
detected_aeq   = f1["ODPT.aeq"][()]
telescope_lamda= f1["ODPT.latitude"][()]
population     = f1["population"][()]
t              = f1["t"][()]
Ekev0          = f1["Ekev0"][()]
precip_id      = f1["precip_id"][()]
precip_lamda   = f1["precip_lamda"][()]
precip_aeq     = f1["precip_aeq"][()]
precip_alpha   = f1["precip_alpha"][()]
precip_time    = f1["precip_time"][()] 
neg_id         = f1["neg_id"][()]
high_id        = f1["high_id"][()]
nan_id         = f1["nan_id"][()]
#noWPI and WPI afterwards read
#saved_id    = f1["saved_id"][()]
#saved_lamda = f1["saved_lamda"][()]
#saved_aeq   = f1["saved_aeq"][()]
#saved_time  = f1["saved_time"][()]
#saved_alpha = f1["saved_alpha"][()]
#saved_ppar  = f1["saved_ppar"][()]
#saved_pper  = f1["saved_pper"][()]
f1.close()

f2 = h5py.File("h5files/20000p_both.h5","r")
#print("Keys: %s" % f2.keys())
detected_lamda_both = f2["ODPT.lamda"][()]
detected_time_both  = f2["ODPT.time"][()]
detected_id_both    = f2["ODPT.id"][()]
detected_alpha_both = f2["ODPT.alpha"][()]
detected_aeq_both   = f2["ODPT.aeq"][()]
telescope_lamda_both= f2["ODPT.latitude"][()]
population_both     = f2["population"][()]
t_both              = f2["t"][()]
Ekev0_both          = f2["Ekev0"][()]
precip_id_both      = f2["precip_id"][()]
precip_lamda_both   = f2["precip_lamda"][()]
precip_aeq_both     = f2["precip_aeq"][()]
precip_alpha_both   = f2["precip_alpha"][()]
precip_time_both    = f2["precip_time"][()] 
neg_id_both         = f2["neg_id"][()]
high_id_both        = f2["high_id"][()]
nan_id_both         = f2["nan_id"][()]
#noWPI and WPI afterwards read
#savedwpi_id    = f2["saved_id"][()]
#savedwpi_lamda = f2["saved_lamda"][()]
#savedwpi_aeq   = f2["saved_aeq"][()]
#savedwpi_time  = f2["saved_time"][()]
#savedwpi_alpha = f2["saved_alpha"][()]
#savedwpi_ppar  = f2["saved_ppar"][()]
#savedwpi_pper  = f2["saved_pper"][()]
f2.close()

############################# TELESCOPE SPECIFICATION && VARIABLES #######################################
time_bin  = 2                #seconds to distinquish events(time resolution)
timesteps = math.ceil(t / time_bin) # t stops at the last timestep (e.g 14.9)
view = 180
sector_range = 2 #P.A bins #1deg
sectors = int(view/sector_range)
##########################################################################################################
##########################################################################################################
###################################### POST PROCESSING - PLOTS ###########################################
##########################################################################################################
##########################################################################################################
print()
print("Population     :",population)
print("Simulation time:",t,"seconds")
print("Time bins      :",time_bin,"   seconds") #Sampling frequency should comply with Nyquist
print("Timesteps      :",timesteps) 
print("P.A bins       :",sector_range,"degrees")
print("Sectors        :",sectors)
print()
print("|noWPI   WPI|")
print(" ",len(precip_id),"     ",len(precip_id_both)," particles escaped")
print(" ",len(neg_id),"     ",len(neg_id_both),"particles developed negative P.A")
print(" ",len(high_id),"     ",len(high_id_both),"  particles developed high P.A")
print(" ",len(nan_id),"     ",len(nan_id_both),"  particles developed nan P.A") 
print("These particles will be excluded from the bining population losing:", ((len(neg_id)+len(neg_id_both)+len(nan_id)+len(nan_id_both)+len(high_id)+len(high_id_both))/population)*100,"% of the population\n")


########################################### FONTS AND COLORS #############################################
font = {'family': 'serif',
        'color':  'blue',
        'weight': 'bold',
        'size': 12}
colors = []
for i in range(max(timesteps,sectors)):      #colors to seperate timesteps or sectors.
    colors.append('#%06X' % randint(0, 0xFFFFFF))
######################################## PLOT SAVED PARTICLE #######################################
"""$fig,ax = plt.subplots(2,2)$ax[0,0].scatter(saved_time,saved_lamda*R2D, s=1)$ax[0,0].scatter(savedwpi_time,savedwpi_lamda*R2D, s=1, alpha=0.01)$ax[0,0].set(xlabel="Time", ylabel="lamda")$ax[0,1].scatter(saved_lamda*R2D,saved_alpha*R2D, s=1)$ax[0,1].scatter(savedwpi_lamda*R2D,savedwpi_alpha*R2D, s=1, alpha=0.01)$ax[0,1].set(xlabel="lamda", ylabel="alpha")$ax[1,0].scatter(saved_alpha*R2D,saved_ppar, s=1)$ax[1,0].scatter(savedwpi_alpha*R2D,savedwpi_ppar, s=1, alpha=0.01)$ax[1,0].set(xlabel="alpha", ylabel="ppar")$ax[1,1].scatter(saved_alpha*R2D,saved_pper, s=1)$ax[1,1].scatter(savedwpi_alpha*R2D,savedwpi_pper, s=1, alpha=0.01)$ax[1,1].set(xlabel="alpha", ylabel="pper")$fig.savefig("SingleParticle.png",dpi=200)"""
################################### CROSSING PARTICLES LAMDA-TIME PLOT ####################################

#noWPI
fig, ax = plt.subplots()
ax.scatter(detected_time, detected_lamda*R2D, c = detected_id, s=0.3, cmap="viridis")
ax.grid(alpha=0.8)
ax.set(xlabel="time(s)", ylabel="latitude(deg)")
plt.title("$Population$: " +str(population)+"\nNorthward particles are captured below the satellite.\nSouthward particles are captured above the satellite",size="medium")
plt.annotate("SATELLITE",xy=(t/2,telescope_lamda+0.0002),color="blue",weight="semibold")
ax.ticklabel_format(useOffset=False)    #disable e notation.
ax.axhline(y = telescope_lamda ,color="b", linestyle="dashed")
plt.savefig("simulation_MM/Crossing_particles.png", dpi=100)
#BOTH
fig, ax = plt.subplots()
ax.scatter(detected_time_both, detected_lamda_both*R2D, c = detected_id_both, s=0.3, cmap="viridis")
ax.grid(alpha=0.8)
ax.set(xlabel="time(s)", ylabel="latitude(deg)")
plt.title("$Population$: " +str(population_both)+"\nNorthward particles are captured below the satellite.\nSouthward particles are captured above the satellite",size="medium")
plt.annotate("SATELLITE",xy=(t_both/2,telescope_lamda_both+0.0002),color="blue",weight="semibold")
ax.ticklabel_format(useOffset=False)    #disable e notation.
ax.axhline(y = telescope_lamda_both ,color="b", linestyle="dashed")
plt.savefig("simulation_MM/Crossing_particles_both.png", dpi=100)

############################################## BINNING ####################################################
#"""
#noWPI
sctr_flux = [ [0 for i in range(timesteps)] for j in range(sectors) ]   #sctr_flux[sectors][timesteps]
sum_flux  = [  0 for i in range(timesteps)]

###BINNING WITH LOCAL PITCH ANGLE: detected_alpha not detected_aeq
#Formula for aeq is not valid. Do the binning with local P.A since satellite is @0deg aeq~=alpha.
for time,pa,id in zip(detected_time,detected_alpha,detected_id): #Iterate in both array elements. No sorting required.
    if (pa<0):
        print("Negative pa")
        sys.exit(1)
    timestep = math.floor(time/time_bin)
    sector   = math.floor(pa*R2D/sector_range)
    if(sector>=sectors or timestep>=timesteps): #Uneeded? Happens when NAN particles are kept in the simulation(jumping steps) or when initializing close to 0 and 180 aeq. WHY?
        print("noWPI Particle",id,"at time",time,"needs sector",sector,"in timestep",timestep)
        #sys.exit(1) 
        continue 
    if( (id in neg_id) or (id in high_id) or (id in nan_id) or (id in neg_id_both) or (id in high_id_both) or (id in nan_id_both)):
        continue #don't include these particles
    if(pa*R2D==180):
        sector = sectors-1 #to include p.a 180 in the last sector. Is this needed?
    sctr_flux[sector][timestep] += 1              #Number of detected particles in this sector-timestep.
    sum_flux[timestep] += 1

#BOTH
sctr_flux_both = [ [0 for i in range(timesteps)] for j in range(sectors) ]   #sctr_flux[sectors][timesteps]
sum_flux_both  = [  0 for i in range(timesteps)]
for time,pa,id in zip(detected_time_both,detected_alpha_both,detected_id_both):
    if (pa<0):
        print("Negative pa")
        sys.exit(1)
    timestep = math.floor(time/time_bin)
    sector   = math.floor(pa*R2D/sector_range)
    #print("WPI Particle",id,"at time",time,"needs sector",sector,"in timestep",timestep)
    if(sector>=sectors or timestep>=timesteps): #Uneeded? Happens when NAN particles are kept in the simulation(jumping steps) or when initializing close to 0 and 180 aeq. WHY?
        print("WPI Particle",id,"with pa",pa*R2D,"at time",time,"needs sector",sector,"in timestep",timestep)
        #sys.exit(1)
        continue 
    if( (id in neg_id) or (id in high_id) or (id in nan_id) or (id in neg_id_both) or (id in high_id_both) or (id in nan_id_both)):
        continue #don't include these particles
    if(pa*R2D==180):
        sector = sectors-1 #to include p.a 180 in the last sector. Is this needed?
    sctr_flux_both[sector][timestep] += 1              #Number of detected particles in this sector-timestep.
    sum_flux_both[timestep] += 1

#PRECIPITATING
#These particles escape with alpha close to 90(?), pitch angle binning to compare them with the particles that are crossing the equator
sctr_flux_precip = [ [0 for i in range(timesteps)] for j in range(sectors) ]   #sctr_flux[sectors][timesteps]
sum_flux_precip = [0 for i in range(timesteps)]
for time,pa in zip(precip_time_both,precip_aeq_both):
    if (pa<0):
        print("Negative pa")
        sys.exit(1)
    timestep = math.floor(time/time_bin)
    sector   = math.floor(pa*R2D/sector_range)
    if(sector>=sectors or timestep>=timesteps): #Uneeded? Happens when NAN particles are kept in the simulation(jumping steps) or when initializing close to 0 and 180 aeq. WHY?
        print("Precipitating Particle at time",time,"needs sector",sector,"in timestep",timestep)
        #sys.exit(1) 
        continue 
    if(pa*R2D==180):
        sector = sectors-1 #to include p.a 180 in the last sector. Is this needed?
    sctr_flux_precip[sector][timestep] += 1              #Number of detected particles in this sector-timestep.
    sum_flux_precip[timestep] += 1
######################################### PARTICLE SUM - 360 PLOT ###########################################
#"""
fig, ax = plt.subplots()
plt.title("Detected particle sum in all look_dirs for "+str(t)+" seconds, in "+str(timesteps)+" timesteps\n Satellite @"+str(telescope_lamda)+" deg")
ax.set(xlabel="Time(s), in time_bins of "+str(time_bin)+"(s)", ylabel="Total Flux")
#ax.set_xticks(ticks=np.arange(0,t+time_bin,time_bin),minor=True) #ticks for time_bin seperation
ax.xaxis.grid(True, which='both')
for timestep in range(0,timesteps):
    if sum_flux[timestep]!=0:
        ax.scatter(timestep*time_bin+time_bin/2,sum_flux[timestep],s=10, c="black",label="noWPI")
    if sum_flux_both[timestep]!=0:  
        ax.scatter(timestep*time_bin+time_bin/2,sum_flux_both[timestep],s=10, c="red", label="noWPI && WPI")  
    if timestep==0:#plot legend once
        ax.legend()
plt.savefig("simulation_MM/Particle_sum_bins"+str(time_bin)+"s.png", dpi=100)
#"""
##################################### WPI-NOWPI DIFF FOR HISTOGRAM #########################################
#"""
moved  = [ [0 for i in range(sectors)] for j in range(timesteps) ]   #sctr_flux[timesteps][sectors]
precip = [ [0 for i in range(sectors)] for j in range(timesteps) ]   
lost   = [ 0 for i in range(sectors)]   

for sector in range(0,sectors):  
    for timestep in range(0,timesteps): 
        precip[timestep][sector] = sctr_flux_precip[sector][timestep] #particles that got lost in this timestep and sector.
        lost[sector]             = lost[sector] + precip[timestep][sector]        #lost particles until this timestep, for this sector.
        moved[timestep][sector]  = abs(sctr_flux[sector][timestep] - sctr_flux_both[sector][timestep]) - lost[sector] #find difference between simulations and remove lost particles
 #"""       
###################################### (FLUX-P.A)*TIMESTEPS MOVIE ##########################################
#"""
fig,ax = plt.subplots()
FFMpegWriter = manimation.writers["ffmpeg"]
metadata2 = dict(title="P.A binning",comment="P.A bins of"+str(sector_range)+"degrees")
fps = 1
writer = FFMpegWriter(fps=fps, metadata = metadata2)
print("Generating P.A binning mp4 file...\nDuration of mp4 file will be:",(timesteps*fps), "seconds")
with writer.saving(fig, "simulation_MM/Bins_"+str(sector_range)+"deg_"+str(time_bin)+"s.mp4", 100):
    for timestep in range(0,timesteps):           
            
        for sector in range(0,sectors): 

            if(sctr_flux[sector][timestep]!=0): #plot only non zero dots
                ax.scatter(sector+0.5, sctr_flux[sector][timestep],c="black",s=1) 
            if(sctr_flux_both[sector][timestep]!=0 and sctr_flux_both[sector][timestep]!=sctr_flux[sector][timestep]): #plot dot only when it differs with noWPI
                ax.scatter(sector+0.5, sctr_flux_both[sector][timestep],c="red",s=1) 
             
        #P.A change and precipitation histograms
        ax.bar(np.arange(0.5,sectors+0.5),moved[timestep],  width=0.3, color='blue',label="moved particles")               #plot difference with bars
        ax.bar(np.arange(0.5,sectors+0.5),precip[timestep], width=0.3, bottom=moved[timestep], color='orange', label="precipitated particles")#plot difference precipitated
        
        
        ax.set_yscale("log")
        #Sector red ticks
        ax.set_xticks(ticks=np.arange(0,sectors,10)) 
        ax.set_xticklabels(labels=np.arange(0,sectors,10),color="red",size="small")
        #P.A black ticks
        #ax.set_xticks(ticks=np.arange(0,view,10),minor=True) 
        #ax.set_xticklabels(labels=np.arange(0,view,10),minor=True,size="small") 
        #Limits
        ax.set_ylim(0, max(np.amax(sctr_flux), np.amax(sctr_flux_both)) )
        ax.set_xlim(0,sectors)

        ax.set(ylabel="count")
        ax.set(xlabel="P.A bins(black): "+str(sector_range)+" deg   Sectors(red): "+str(sectors))
        ax.set_title("Equatorial P.A dstr, $time: "+str("{:.1f}".format(timestep*time_bin))+"s$", loc="left", size="small",color="blue",x=-0.15)
        ax.xaxis.grid(True, which='major')
        ax.xaxis.grid(True, which='minor')
        ax.legend(loc='upper right')
        writer.grab_frame()
        ax.clear() #clear data 
#"""
