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

############################################# READ HDF5 ###################################################
#noWPI read
f1 = h5py.File("h5files/1p_no_wpi.h5","r")
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
lamda00        = f1["lamda00"][()] #States when noWPI stops
ppar00         = f1["ppar00"][()]
pper00         = f1["pper00"][()]
alpha00        = f1["alpha00"][()]
aeq00          = f1["aeq00"][()]
eta00          = f1["eta00"][()]
time00         = f1["time00"][()]
precip_id      = f1["precip_id"][()]
precip_lamda   = f1["precip_lamda"][()]
precip_alpha   = f1["precip_alpha"][()]
precip_aeq     = f1["precip_aeq"][()]
precip_time    = f1["precip_time"][()] 
saved_id    = f1["saved_id"][()]
saved_lamda = f1["saved_lamda"][()]
saved_alpha = f1["saved_alpha"][()]
saved_aeq   = f1["saved_aeq"][()]
saved_ppar  = f1["saved_ppar"][()]
saved_pper  = f1["saved_pper"][()]
saved_time  = f1["saved_time"][()]
f1.close()

#noWPI and WPI afterwards read
f2 = h5py.File("h5files/1p_both.h5","r")
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
lamda00_both        = f2["lamda00"][()] #States when noWPI stops
ppar00_both         = f2["ppar00"][()]
pper00_both         = f2["pper00"][()]
alpha00_both        = f2["alpha00"][()]
aeq00_both          = f2["aeq00"][()]
eta00_both          = f2["eta00"][()]
time00_both         = f2["time00"][()]
precip_id_both      = f2["precip_id"][()]
precip_lamda_both   = f2["precip_lamda"][()]
precip_alpha_both   = f2["precip_alpha"][()]
precip_aeq_both     = f2["precip_aeq"][()]
precip_time_both    = f2["precip_time"][()] 
savedwpi_id    = f2["saved_id"][()]
savedwpi_lamda = f2["saved_lamda"][()]
savedwpi_alpha = f2["saved_alpha"][()]
savedwpi_aeq   = f2["saved_aeq"][()]
savedwpi_ppar  = f2["saved_ppar"][()]
savedwpi_pper  = f2["saved_pper"][()]
savedwpi_time  = f2["saved_time"][()]
f2.close()

#Distribution read
f3 = h5py.File("h5files/1p_testAEQ_testLAMDA.h5","r")
aeq0         = f3["aeq"][()]
lamda0       = f3["lat"][()]
aeq0_bins    = f3["aeq0_bins"][()]
f3.close()
#Save initials to CSV file
header = ['id', 'aeq0', 'lamda0', 'aeq00 after nowpi', 'lamda00 after nowpi']
data = []
p=0
for a0,l0,a00,l00 in zip(aeq0, lamda0, aeq00_both, lamda00_both):
    data.extend([[p,a0*R2D,l0*R2D,a00*R2D,l00*R2D]])
    p=p+1 #id is the element's location in the list since the particles were saved like that
with open("dstr_before_and_after_WPI.csv", "w") as f:
    writer = csv.writer(f)
    writer.writerow(header)
    writer.writerows(data)
##########################################################################################################
##########################################################################################################
###################################### POST PROCESSING - PLOTS ###########################################
##########################################################################################################
##########################################################################################################




############################# TELESCOPE SPECIFICATION && VARIABLES #######################################
time_bin  = 3                 #seconds to distinquish events(time resolution)
timesteps = int (t / time_bin)
view = 180 
sector_range = 10 #P.A bins #1deg
sectors = int(view/sector_range)
########################################### FONTS AND COLORS #############################################
font = {'family': 'serif',
        'color':  'blue',
        'weight': 'bold',
        'size': 12}
colors = []
for i in range(max(timesteps,sectors)):      #colors to seperate timesteps or sectors.
    colors.append('#%06X' % randint(0, 0xFFFFFF))


######################################## PLOT SAVED PARTICLE #######################################
#"""
fig,ax = plt.subplots(2,2)
ax[0,0].scatter(saved_time,saved_lamda*R2D, s=1)
ax[0,0].scatter(savedwpi_time,savedwpi_lamda*R2D, s=1, alpha=0.01)
ax[0,0].set(xlabel="Time", ylabel="lamda")

ax[0,1].scatter(saved_lamda*R2D,saved_alpha*R2D, s=1)
ax[0,1].scatter(savedwpi_lamda*R2D,savedwpi_alpha*R2D, s=1, alpha=0.01)
ax[0,1].set(xlabel="lamda", ylabel="alpha")

ax[1,0].scatter(saved_alpha*R2D,saved_ppar, s=1)
ax[1,0].scatter(savedwpi_alpha*R2D,savedwpi_ppar, s=1, alpha=0.01)
ax[1,0].set(xlabel="alpha", ylabel="ppar")

ax[1,1].scatter(saved_alpha*R2D,saved_pper, s=1)
ax[1,1].scatter(savedwpi_alpha*R2D,savedwpi_pper, s=1, alpha=0.01)
ax[1,1].set(xlabel="alpha", ylabel="pper")
#fig.savefig("SingleParticle.png",dpi=200)

#"""
######################################## PLOT INITIAL DISTRIBUTION #######################################
#"""
fig, ax = plt.subplots()
for sec in range(0,sectors):
    ax.scatter(sec,aeq0_bins[sec],s=2,alpha=1)
ax.grid(alpha=.3)
ax.set(xlabel="Sectors",xlim=(0,sectors-1),xticks=np.arange(0,sectors),ylabel="dN",title="Aeq0 distribution, sector range "+str(sector_range)+" degrees")
ax.set_yscale("log")
fig.savefig("simulation_MM/10e6_normals.png",dpi=200)

fig, ax = plt.subplots()
ax.scatter(lamda0*R2D,aeq0*R2D,s=0.5,alpha=0.1)
ax.grid(alpha=.3)
ax.set(xlabel="Latitude(deg)",ylabel="Equatorial P.A",title="Initial lat-aeq of simulated particles",ylim=(1,179),xlim=(-90,90),xticks=np.linspace(-90,90,5))
ax.axhline(y = 90, color ="b", linestyle="dashed")
fig.savefig("simulation_MM/10e6_normals_2.png",dpi=200)
#"""
################################### CROSSING PARTICLES LAMDA-TIME PLOT ####################################
#"""
#noWPI
fig, ax = plt.subplots()
ax.scatter(detected_time, detected_lamda*R2D, c = detected_id, s=0.3, cmap="viridis")
ax.grid(alpha=0.8)
ax.set(xlabel="time(s)", ylabel="latitude(deg)")
plt.title("$Population$: " +str(population)+"]\nNorthward particles are captured below the satellite.\nSouthward particles are captured above the satellite",size="medium")
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
#"""
############################################## BINNING ####################################################
#"""
#noWPI
sctr_flux = [ [0 for i in range(timesteps)] for j in range(sectors) ]   #sctr_flux[sectors][timesteps]
sum_flux  = [  0 for i in range(timesteps)]
for time,pa in zip(detected_time,detected_aeq): #Iterate in both array elements. No sorting required.
    timestep = math.floor(time/time_bin)
    sector   = math.floor(pa*R2D/sector_range)
    if(sector>=sectors or timestep>=timesteps): #Uneeded? Happens when NAN particles are kept in the simulation(jumping steps) or when initializing close to 0 and 180 aeq. WHY?
        print("Particle",id,"with pa",pa*R2D,"needs sector",sector,"in timestep",timestep)
        continue 
    if(pa*R2D==180):
        sector = sectors-1 #to include p.a 180 in the last sector. Is this needed?
    sctr_flux[sector][timestep] += 1              #Number of detected particles in this sector-timestep.
    sum_flux[timestep] += 1

#BOTH
sctr_flux_both = [ [0 for i in range(timesteps)] for j in range(sectors) ]   #sctr_flux[sectors][timesteps]
sum_flux_both  = [  0 for i in range(timesteps)]
for time,pa,id in zip(detected_time_both,detected_aeq_both,detected_id):
    timestep = math.floor(time/time_bin)
    sector   = math.floor(pa*R2D/sector_range)
    if(sector>=sectors or timestep>=timesteps): #Uneeded? Happens when NAN particles are kept in the simulation(jumping steps) or when initializing close to 0 and 180 aeq. WHY?
        print("Particle",id,"with pa",pa*R2D,"needs sector",sector,"in timestep",timestep)
        continue 
    if(pa*R2D==180):
        sector = sectors-1 #to include p.a 180 in the last sector. Is this needed?
    sctr_flux_both[sector][timestep] += 1              #Number of detected particles in this sector-timestep.
    sum_flux_both[timestep] += 1

#PRECIPITATING
#These particles escape with alpha close to 90(?), binning with equatorial pitch angle to compare them with the particles that are crossing the equator--> equatorial P.A
sctr_flux_precip = [ [0 for i in range(timesteps)] for j in range(sectors) ]   #sctr_flux[sectors][timesteps]
sum_flux_precip = [0 for i in range(timesteps)]
for time,pa in zip(precip_time_both,precip_aeq_both):
    timestep = math.floor(time/time_bin)
    sector   = math.floor(pa*R2D/sector_range)
    #if(sector>=sectors or timestep>=timesteps): #Uneeded? Happens when NAN particles are kept in the simulation(jumping steps) or when initializing close to 0 and 180 aeq. WHY?
    #    print("Particle with pa",pa*R2D,"needs sector",sector,"in timestep",timestep)
    #    continue 
    if(pa*R2D==180):
        sector = sectors-1 #to include p.a 180 in the last sector. Is this needed?
    sctr_flux_precip[sector][timestep] += 1              #Number of detected particles in this sector-timestep.
    sum_flux_precip[timestep] += 1
#"""
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
        lost[sector] = lost[sector] + precip[timestep][sector]        #lost particles until this timestep, for this sector.
        moved[timestep][sector] = abs(sctr_flux[sector][timestep] - sctr_flux_both[sector][timestep]) - lost[sector] #find difference between simulations and remove lost particles
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
    for timestep in range(0,timestep):           
            
        for sector in range(0,sectors): 

            if(sctr_flux[sector][timestep]!=0): #plot only non zero dots
                ax.scatter(sector+0.5, sctr_flux[sector][timestep],c="black",s=1) 
            if(sctr_flux_both[sector][timestep]!=sctr_flux[sector][timestep]): #plot dot only when it differs with noWPI
                ax.scatter(sector+0.5, sctr_flux_both[sector][timestep],c="red",s=1) 

        ax.bar(np.arange(0.5,sectors+0.5),moved[timestep], width=0.3,color='blue',label="moved particles")               #plot difference with bars
        ax.bar(np.arange(0.5,sectors+0.5),precip[timestep], width=0.3, bottom=moved[timestep], color='orange', label="precipitated particles")#plot difference precipitated
        ax.set_yscale("log")

        #Sector red ticks
        ax.set_xticks(ticks=np.arange(0.5,sectors,10)) 
        ax.set_xticklabels(labels=np.arange(0,sectors,10),color="red",size="small")
        #P.A black ticks
        #ax.set_xticks(ticks=np.arange(0,sectors,),minor=True) 
        #ax.set_xticklabels(labels=np.arange(0,view,sector_range),minor=True,size="small") 

        ax.set_ylim(0, max(np.amax(sctr_flux)+(0.1*np.amax(sctr_flux)), np.amax(sctr_flux_both)+(0.1*np.amax(sctr_flux_both)) ) )
        ax.set_xlim(0,sectors)
        ax.set(ylabel="count")
        ax.set(xlabel="P.A bins(black): "+str(sector_range)+" deg   Sectors(red): "+str(sectors))
        ax.set_title("Equatorial P.A dstr, $time: "+str("{:.1f}".format(timestep*time_bin))+"s$", loc="left", size="small",color="blue",x=-0.15)
        ax.xaxis.grid(True, which='major')
        ax.xaxis.grid(True, which='minor')
        ax.legend()
        writer.grab_frame()
        ax.clear() #clear data 
#"""