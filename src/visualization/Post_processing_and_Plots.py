import numpy as np 
import sys
import matplotlib.pyplot as plt
import h5py
import math
import matplotlib.animation as manimation
from random import randint
import threading
import os

np.set_printoptions(threshold=sys.maxsize)
D2R=np.pi/180
R2D=1/D2R
############################################# READ HDF5 ###################################################
#noWPI read
f1 = h5py.File("output/files/1p_nowpi.h5","r")
detected_lamda = f1["ODPT.lamda"][()]
detected_time  = f1["ODPT.time"][()]
detected_id    = f1["ODPT.id"][()]
detected_alpha = f1["ODPT.alpha"][()]
detected_aeq   = f1["ODPT.aeq"][()]
telescope_lamda= f1["ODPT.latitude"][()]
population     = f1["population"][()]
t              = f1["t"][()]
precip_id      = f1["precip_id"][()]
precip_lamda   = f1["precip_lamda"][()]
precip_aeq     = f1["precip_aeq"][()]
precip_alpha   = f1["precip_alpha"][()]
precip_time    = f1["precip_time"][()] 
neg_id         = f1["neg_id"][()]
high_id        = f1["high_id"][()]
nan_id         = f1["nan_id"][()]
f1.close()
"""
f2 = h5py.File("output/files/50000p_nowpi_wpi.h5","r")
detected_lamda_both = f2["ODPT.lamda"][()]
detected_time_both  = f2["ODPT.time"][()]
detected_id_both    = f2["ODPT.id"][()]
detected_alpha_both = f2["ODPT.alpha"][()]
detected_aeq_both   = f2["ODPT.aeq"][()]
telescope_lamda_both= f2["ODPT.latitude"][()]
population_both     = f2["population"][()]
t_both              = f2["t"][()]
precip_id_both      = f2["precip_id"][()]
precip_lamda_both   = f2["precip_lamda"][()]
precip_aeq_both     = f2["precip_aeq"][()]
precip_alpha_both   = f2["precip_alpha"][()]
precip_time_both    = f2["precip_time"][()] 
neg_id_both         = f2["neg_id"][()]
high_id_both        = f2["high_id"][()]
nan_id_both         = f2["nan_id"][()]
saved_id            = f2["saved_id"][()]
saved_deta_dt       = f2["saved_deta_dt"][()]
saved_lamda         = f2["saved_lamda"][()]
saved_Ekin          = f2["saved_Ekin"][()]
f2.close()

######################### TELESCOPE SPECIFICATION -- BINNING PARAMETERS ###############################
time_bin  = 2                #seconds to distinquish events(time resolution)
timesteps = math.ceil(t / time_bin) # t stops at the last timestep (e.g 14.9)
view = 180
sector_range = 2 #P.A bins, look directions
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
#Create directory for plots"""
filepath_plots = 'output/plots/'+str(population)+"p_"+str(int(t))+"s"
if (os.path.exists(filepath_plots)):
    print ("The directory %s already exists" % filepath_plots)
else:
    try:
        os.makedirs(filepath_plots)
    except OSError as error:
        print(error)
        print ("Creation of the directory %s failed, " % filepath_plots)
    else:
        print ("Successfully created the directory %s" % filepath_plots)
"""
########################################### FONTS AND COLORS #############################################
font = {'family': 'serif',
        'color':  'blue',
        'weight': 'bold',
        'size': 12}
colors = []
for i in range(max(timesteps,sectors)):      #colors to seperate timesteps or sectors.
    colors.append('#%06X' % randint(0, 0xFFFFFF))"""
################################### CROSSING PARTICLES LAMDA-TIME PLOT ####################################
#noWPI
#"""
fig, ax = plt.subplots()
ax.scatter(detected_time, detected_lamda*R2D, c = detected_id, s=10, cmap="viridis")
ax.grid(alpha=0.8)
ax.set(xlabel="time(s)", ylabel="latitude(deg)")
plt.title("$Population$: " +str(population)+"\nNorthward particles are captured below the satellite.\nSouthward particles are captured above the satellite",size="medium")
plt.annotate("SATELLITE",xy=(t/2,telescope_lamda+0.0002),color="blue",weight="semibold")
ax.ticklabel_format(useOffset=False)    #disable e notation.
ax.axhline(y = telescope_lamda ,color="b", linestyle="dashed")
plt.savefig(filepath_plots+"/Crossing_particles.png", dpi=100)
#BOTH
fig, ax = plt.subplots()
ax.scatter(detected_time_both, detected_lamda_both*R2D, c = detected_id_both, s=0.3, cmap="viridis")
ax.grid(alpha=0.8)
ax.set(xlabel="time(s)", ylabel="latitude(deg)")
plt.title("$Population$: " +str(population_both)+"\nNorthward particles are captured below the satellite.\nSouthward particles are captured above the satellite",size="medium")
plt.annotate("SATELLITE",xy=(t_both/2,telescope_lamda_both+0.0002),color="blue",weight="semibold")
ax.ticklabel_format(useOffset=False)    #disable e notation.
ax.axhline(y = telescope_lamda_both ,color="b", linestyle="dashed")
plt.savefig(filepath_plots+"/Crossing_particles_both.png", dpi=100)
#"""

########################################## MINIMUM DETA_DT ################################################
#For every particle in the population the state for minimum deta_dt is saved. 
fix,ax = plt.subplots()
cc = ax.scatter(np.rad2deg(saved_lamda),saved_deta_dt,s=1, c=saved_Ekin,cmap='jet')
ax.axhline(y = 0 ,color="b", linestyle="dashed")
ax.set(xlabel="latitude[deg]",ylabel="deta_dt")
cbar=plt.colorbar(cc)
cbar.ax.set_title('Ekin[keV]')
plt.savefig(filepath_plots+"/min_deta_dt.png")
########################################## MINIMUM DETA_DT ################################################


############################################## BINNING ####################################################
###BINNING WITH LOCAL PITCH ANGLE: detected_alpha not detected_aeq
#Formula for aeq is not valid. Do the binning with local P.A since satellite is @0deg aeq~=alpha.
#"""
sctr_flux = [ [0 for i in range(timesteps)] for j in range(sectors) ]   #sctr_flux[sectors][timesteps]
sum_flux  = [  0 for i in range(timesteps)]
sctr_flux_both = [ [0 for i in range(timesteps)] for j in range(sectors) ]   
sum_flux_both  = [  0 for i in range(timesteps)]
sctr_flux_precip = [ [0 for i in range(timesteps)] for j in range(sectors) ]   
sum_flux_precip = [0 for i in range(timesteps)]


def thread_noWPI_binning():
#noWPI
    for time,pa,id in zip(detected_time,detected_alpha,detected_id): #Iterate in both array elements. No sorting required.
        timestep = math.floor(time/time_bin)
        sector   = math.floor(pa*R2D/sector_range)
        if (pa<0):
            raise Exception("Negative pa")
        if(sector>=sectors or timestep>=timesteps): #Uneeded? Happens when NAN particles are kept in the simulation(jumping steps) or when initializing close to 0 and 180 aeq. WHY?
            raise Exception("noWPI Particle",id,"at time",time,"needs sector",sector,"in timestep",timestep)
        if( (id in neg_id) or (id in high_id) or (id in nan_id) or (id in neg_id_both) or (id in high_id_both) or (id in nan_id_both)):
            continue #don't include these particles in the binning process
        if(pa*R2D==180):
            sector = sectors-1 #to include particles with p.a==180 in the last sector. Is this needed?
        sctr_flux[sector][timestep] += 1              #Number of detected particles in this sector-timestep.
        sum_flux[timestep] += 1                       #Sum of detected particles in this timestep from all look directions.
    print("Binning noWPI done!")

def thread_both_binning():
#noWPI and then WPI
    for time,pa,id in zip(detected_time_both,detected_alpha_both,detected_id_both):
        timestep = math.floor(time/time_bin)
        sector   = math.floor(pa*R2D/sector_range)
        if (pa<0):
            raise Exception("Negative pa")
        if(sector>=sectors or timestep>=timesteps): 
            raise Exception("WPI Particle",id,"with pa",pa*R2D,"at time",time,"needs sector",sector,"in timestep",timestep)
        if( (id in neg_id) or (id in high_id) or (id in nan_id) or (id in neg_id_both) or (id in high_id_both) or (id in nan_id_both)):
            continue
        if(pa*R2D==180):
            sector = sectors-1 
        sctr_flux_both[sector][timestep] += 1              
        sum_flux_both[timestep] += 1
    print("Binning noWPI+WPI done!")

def thread_precip_binning():
#Prepitating. Just iterating for the latter simulation(noWPI and then WPI) is enough. Particles won't precipitate after one bounce in the noWPI simulation.
#These particles escape with alpha close to 90(?), pitch angle binning to compare them with the particles that are crossing the equator
    for time,pa in zip(precip_time_both,precip_aeq_both):
        timestep = math.floor(time/time_bin)
        sector   = math.floor(pa*R2D/sector_range)
        if (pa<0):
            raise Exception("Negative pa")
        if(sector>=sectors or timestep>=timesteps): 
            raise Exception("Precipitating Particle at time",time,"needs sector",sector,"in timestep",timestep)
        if( (id in neg_id) or (id in high_id) or (id in nan_id) or (id in neg_id_both) or (id in high_id_both) or (id in nan_id_both)):
            continue
        if(pa*R2D==180):
            sector = sectors-1 
        sctr_flux_precip[sector][timestep] += 1           
        sum_flux_precip[timestep] += 1
    print("Binning Precipitating done!")

#Create thread classes
th1=threading.Thread(target=thread_noWPI_binning)
th2=threading.Thread(target=thread_both_binning)
th3=threading.Thread(target=thread_precip_binning)
th1.start()
th2.start()
th3.start()

th1.join()
th2.join()
th3.join()
#main waits until all threads are done

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
plt.savefig(filepath_plots+"/Particle_sum_bins"+str(time_bin)+"s.png", dpi=100)
#"""
##################################### WPI-NOWPI DIFF FOR HISTOGRAM #########################################
#"""
moved  = [ [0 for i in range(sectors)] for j in range(timesteps) ]   #sctr_flux[timesteps][sectors]
precip = [ [0 for i in range(sectors)] for j in range(timesteps) ]   
lost   = [ 0 for i in range(sectors)]   

for sector in range(0,sectors):  
    for timestep in range(0,timesteps): 
        precip[timestep][sector] = sctr_flux_precip[sector][timestep] #particles lost in this timestep and sector.
        lost[sector] += precip[timestep][sector]                      #particles lost until this timestep,for this sector.
        moved[timestep][sector]  = abs(sctr_flux[sector][timestep] - sctr_flux_both[sector][timestep]) - lost[sector] #find difference between simulations and remove lost particles until this timestep for this sector
#"""       
###################################### (FLUX-P.A)*TIMESTEPS MOVIE ##########################################
#"""
fig,ax = plt.subplots()
FFMpegWriter = manimation.writers["ffmpeg"]
metadata2 = dict(title="P.A binning",comment="P.A bins of"+str(sector_range)+"degrees")
fps = 1
writer = FFMpegWriter(fps=fps, metadata = metadata2)
print("Generating P.A binning mp4 file...\nDuration of mp4 file will be:",(timesteps*fps), "seconds")
with writer.saving(fig, filepath_plots+"/Bins_"+str(sector_range)+"deg_"+str(time_bin)+"s.mp4", 100):
    for timestep in range(0,timesteps):           
            
        for sector in range(0,sectors): 

            if(sctr_flux[sector][timestep]!=0): #plot only non zero dots
                ax.scatter(sector+0.5, sctr_flux[sector][timestep],c="black",s=1) 
            if(sctr_flux_both[sector][timestep]!=0 and sctr_flux_both[sector][timestep]!=sctr_flux[sector][timestep]): #plot dot only when it differs with noWPI
                ax.scatter(sector+0.5, sctr_flux_both[sector][timestep],c="red",s=1) 
             
        #P.A change and precipitation histograms
        ax.bar(np.arange(0.5,sectors+0.5),moved[timestep],  width=0.3, color='blue',label="moved")               #plot difference with bars
        ax.bar(np.arange(0.5,sectors+0.5),precip[timestep], width=0.3, bottom=moved[timestep], color='orange', label="precip")#plot difference precipitated
        
        
        ax.set_yscale("log")
        ax.set(ylabel="count")
        ax.set(xlabel="P.A bins(black): "+str(sector_range)+" deg   Sectors(red): "+str(sectors))
        ax.set_title("Equatorial P.A dstr, $time: "+str("{:.1f}".format(timestep*time_bin))+"s$", loc="left", size="small",color="blue",x=-0.15)
        ax.set_ylim(1, 100000 )
        ax.set_xlim(0,sectors)
        ax.legend(loc='upper right')

        #Sector red ticks
        ax.set_xticks(ticks=np.arange(0,sectors,10)) 
        ax.set_xticklabels(labels=np.arange(0,sectors,10),color="red",size="small")
        ax.xaxis.grid(True, which='major')
        #P.A black ticks
        #ax.set_xticks(ticks=np.arange(0,view,10),minor=True) 
        #ax.set_xticklabels(labels=np.arange(0,view,10),minor=True,size="small") 
        #ax.xaxis.grid(True, which='minor')

        writer.grab_frame()
        ax.clear() #clear data 
#"""
