import numpy as np 
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


############################################# READ DATA ###################################################
f = h5py.File("h5files/1000p_5s_aeq_time.h5","r")
#print("Keys: %s" % f.keys())
#detected_lamda = f["ODPT.lamda"][()]
#detected_time  = f["ODPT.time"][()]
#detected_id    = f["ODPT.id"][()]
#detected_aeq   = f["ODPT.aeq"][()]
#detected_alpha = f["ODPT.alpha"][()]
#telescope_lamda= f["ODPT.latitude"][()]
#population     = f["population"][()]
#t              = f["t"][()]
#lamda0         = f["lamda0"][()]

aeq            = f["aeq_plot"][()]
alpha          = f["alpha_plot"][()]
time           = f["time_plot"][()]
deta_dt        = f["deta_dt"][()]

f.close();
detected_pop   = len(detected_lamda)

##########################################################################################################
###################################### POST PROCESSING - PLOTS ###########################################
##########################################################################################################

 
def leastFrequent(arr, n):
 
    Hash = dict()
    for i in range(n):
        if arr[i] in Hash.keys():
            Hash[arr[i]] += 1
        else:
            Hash[arr[i]] = 1
 
    min_count = n + 1
    res = -1
    for i in Hash:
        if (min_count >= Hash[i]):
            res = i
            min_count = Hash[i]
         
    return res


least = leastFrequent(detected_id,len(detected_id))


#PLOT PA-TIME
fig, ax = plt.subplots()    
ax.plot(time[least,:-1],alpha[least,:-1])
ax.grid(alpha=.3)
ax.set_xlabel('time (sec)')
ax.set_ylabel('$ a_{eq}$ (deg)')
fig.savefig('a_eq.png' , dpi=200, facecolor='white',transparent=False)

#PLOT DETA_DT-TIME
fig, ax = plt.subplots()    
print(time[least],aeq[least],deta_dt[least])
ax.plot(time[least,:-1],deta_dt[least,:-1])
ax.grid(alpha=.3)
ax.set_xlabel('time (sec)')
ax.set_ylabel('deta_dt')
fig.savefig('deta_dt.png' , dpi=200, facecolor='white',transparent=False)



"""
####################################### TELESCOPE SPECIFICATION ##########################################
time_bin  = 0.1                 #seconds to distinquish events(time resolution)
timesteps = int (t / time_bin)
view = 180 
sector_range = 15
sectors = int(view/sector_range)
######################################### FONT DICT - COLORS  ############################################
font = {'family': 'serif',
        'color':  'blue',
        'weight': 'medium',
        'size': 10,
        }
colors = []
for i in range(max(timesteps,sectors)):      #colors to seperate timesteps or sectors.
    colors.append('#%06X' % randint(0, 0xFFFFFF))
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


################################### CROSSING PARTICLES LAMDA-TIME PLOT ####################################
fig, ax = plt.subplots()


ax.scatter(detected_time, detected_lamda*R2D, c = detected_id, s=1)

ax.grid(alpha=0.8)

ax.set(xlabel="time(s)", ylabel="latitude(deg)")
plt.title("Detected particles in time. [Population: " +str(population)+ ", lamda0: " +str(round(lamda0*R2D))+ "]\n"+"Northward particles are captured below the satellite.\nSouthward particles are captured above the satellite.")
plt.annotate("SATELLITE",xy=(t/2,telescope_lamda+0.0002),color="blue",weight="semibold")
ax.ticklabel_format(useOffset=False)    #disable e notation.
ax.axhline(y = telescope_lamda ,color="b", linestyle="dashed")

plt.savefig("simulation_MM/detected/Crossing_particles.png", dpi=100)

######################################### PARTICLE SUM - 360 PLOT ###########################################
fig, ax = plt.subplots()
plt.title("Detected particle sum in all look_dirs for "+str(t)+" seconds, in "+str(timesteps)+" timesteps\n Satellite @"+str(telescope_lamda)+" deg")
ax.set(xlabel="Time bins-Timesteps in red", ylabel="sum", xticks=np.arange(0,t+time_bin,time_bin))
#deactivate time_bin ticks. no space for t=5s
ax.set_xticklabels(labels=[])
ax.set_xticks(ticks=np.arange(time_bin/2,t,t*time_bin),minor=True)
ax.set_xticklabels(labels=np.arange(0,timesteps,t).astype(int),color="red",minor=True)
ax.grid(alpha=0.8,axis="x")
for timestep in range(0,timesteps):
    if sum_flux[timestep]!=0:
        ax.scatter(timestep*0.1+0.1/2,sum_flux[timestep])

plt.savefig("simulation_MM/detected/Particle_sum.png", dpi=100)
####################################### (FLUX-TIME)*SECTORS  MOVIE ##########################################
fig,ax = plt.subplots()
FFMpegWriter = manimation.writers["ffmpeg"]
metadata = dict(title="Time binning",comment="Time bins of 0.1s")
fps = 1
writer = FFMpegWriter(fps=fps, metadata = metadata)
print("Generating Time binning mp4 file...\nDuration of mp4 file will be:",(sectors*fps), "seconds")


with writer.saving(fig, "simulation_MM/detected/Time_binning.mp4", 100):

    for sector in range(0,sectors):

        for timestep in range(0,timesteps):           
                
            if sctr_flux[sector][timestep] !=0: #to show only non zero dots.
                ax.scatter(timestep*time_bin+(time_bin/2),sctr_flux[sector][timestep],c=colors[timestep])
       
        ax.yaxis.set_major_locator(MaxNLocator(integer=True,nbins=5))
        #ax.set_xticks(ticks=np.arange(0,t+time_bin,time_bin)) 
        ax.set_xticks(ticks=np.arange(0,t+time_bin,time_bin))
        ax.set_xticklabels(labels=[]) #time_bin ticks deactivate for 2+ seconds
        ax.set_xticks(ticks=np.arange(time_bin/2,t,t*time_bin),minor=True) #timestep ticks
        ax.set_xticklabels(labels=np.arange(0,timesteps,t).astype(int),minor=True,c="red")  
        ax.grid(alpha=0.8, axis="x" )
        ax.set_xlabel(xlabel="Time bins of "+str(time_bin)+"s-Timesteps in red")
        ax.set_ylabel(ylabel="Flux")
        ax.set_title("Particle Flux in Sector "+str(sector)+" - Look direction: ("+str(180+sector*sector_range)+"\N{DEGREE SIGN},"+str(180+sector*sector_range+sector_range)+"\N{DEGREE SIGN}) Simulation time: "+str(t)+"s\n                                         -Detectable P.A: [ "+str(sector*sector_range)+"\N{DEGREE SIGN},"+str(sector*sector_range+sector_range)+"\N{DEGREE SIGN}) Satellite @"+str(telescope_lamda)+" deg",loc="left",fontdict=font)              


        writer.grab_frame()
        ax.clear() #clear data 

"""



















"""
###################################### (FLUX-P.A)*TIMESTEPS  MOVIE ##########################################
fig,ax = plt.subplots()
metadata2 = dict(title="P.A binning",comment="P.A bins of 15deg")
fps = 1
writer = FFMpegWriter(fps=fps, metadata = metadata2)
print("Generating P.A binning mp4 file...\nDuration of mp4 file will be:",(timesteps*fps), "seconds")


with writer.saving(fig, "simulation_MM/PA_binning.mp4", 100):

    for timestep in range(0,timesteps):           
            
        for sector in range(13,sectors): #skips 0 sector -7.5 7.5

            ax.set_xlabel("P.A bins") 
            ax.set_ylabel("count") 
            ax.set_title("P.A dstr in Timestep "+str(timestep), loc="left", fontdict=font,x=-0.015)
                        
            if sctr_flux[sector][timestep]!=0: #to show only non zero dots.
                ax.scatter(sector, sctr_flux[sector][timestep],c=colors[sector])  #y axis? flux? p.a can be very close which makes dots merge
                ax.set_xticks(ticks=np.arange(13,sectors)) #pa_bin ticks
                ax.set_xticklabels(labels=np.arange(13,sectors),color="red")
                ax.set_xticks(ticks=np.arange(12.5,sectors,1),minor=True) 
                ax.set_xticklabels(labels=np.arange(7.5,180,15),minor=True,size="small") 
            ax.grid(alpha=0.8, axis="x" )


        writer.grab_frame()
        ax.clear() #clear data 



"""








