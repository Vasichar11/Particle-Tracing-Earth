import numpy as np 
import matplotlib.pyplot as plt
from matplotlib import cm
import h5py
import math
import matplotlib.animation as manimation


D2R=np.pi/180
R2D=1/D2R

################################    READ DATA FOR ALL PARTICLES   #####################################
f1 = h5py.File("particles.h5","r")
#print("Keys: %s" % f1.keys())
lamda    = f1["lamda 2D"][()]
aeq      = f1["aeq 2D"][()]
aeq2     = f1["aeq2 2D"][()]
eta      = f1["eta 2D"][()]
time_sim = f1["simulation time"][()]
upar     = f1["upar"][()]
f1.close();
lamda0 = lamda[:,0]
eta0   = eta[:,0]
aeq0   = aeq[:,0]
eta0_deg = eta0*R2D
aeq0_deg = aeq0*R2D
lamda0_deg = lamda0*R2D
population = len(eta0)
################################# READ DATA FOR DETECTED PARTICLES ####################################
f2 = h5py.File("detected.h5","r")
#print("Keys: %s" % f2.keys())
detected_lamda = f2["detected_lamda 1D"][()]
detected_time  = f2["detected_time 1D"][()]
detected_id    = f2["detected_id 1D"][()]
detected_aeq   = f2["detected_aeq 1D"][()]
detected_alpha = f2["detected_alpha 1D"][()]
telescope_lamda= f2["telescope_lamda scalar"][()]
f2.close();
detected_pop   = len(detected_lamda)
################################# SORT PARTICLES IN TEMPORAL ORDER #####################################
r=0
while r < len(detected_time) - 1:
    
    if(detected_time[r+1]<detected_time[r]):
        
        detected_time[r], detected_time[r+1]  = detected_time[r+1], detected_time[r]   #swap
        detected_lamda[r],detected_lamda[r+1] = detected_lamda[r+1],detected_lamda[r]  
        detected_alpha[r],detected_alpha[r+1] = detected_alpha[r+1],detected_alpha[r]  
        detected_aeq[r],detected_aeq[r+1]     = detected_aeq[r+1],detected_aeq[r]  
        detected_id[r],   detected_id[r+1]    = detected_id[r+1],   detected_id[r]         
        r = -1 #next loop r=0 => sort from start.
    
    r = r+1
#print(detected_time)





##########################################################################################################
###################################### POST PROCESSING - PLOTS ###########################################
##########################################################################################################

####################################### TELESCOPE SPECIFICATION ##########################################
t = round(time_sim[0][-1]) #simulation time
time_bin = 0.1                  #seconds to distinquish events(time resolution)
timesteps = int (t / time_bin)
view = 360 
pa_bin = 10                     #degrees to distinquish incoming particle pitch angles               
pa_sectors = int(180/pa_bin)   
sector_range = 15
sectors = int(view/sector_range)
############################################# FONT DICT  ##################################################
font = {'family': 'serif',
        'color':  'blue',
        'weight': 'medium',
        'size': 10,
        }
colors = np.array(["red","green","blue","yellow","pink","black","orange","purple","beige","brown"]) #to seperate counts in timestep

############################################## BINNING #####################################################
sctr_flux = [ [0 for i in range(timesteps)] for j in range(sectors) ]   #sctr_flux[timesteps][sectors]
pa_flux   = [ [[0 for i in range(pa_sectors)] for j in range(timesteps)] for w in range(sectors) ]   #pa_flux[timesteps][sectors][pa_bins]


for time,pa in zip(detected_time,detected_alpha): #Iterate in both array elements. Particles are sorted in temporal order.
    timestep = math.floor(time*10)
    sector   = math.floor((pa*R2D - 7.5 + 180)/sector_range) + 1 #anagwgi p.a ston sector, sigoura?
    sctr_flux[sector][timestep] += 1                #Number of detected particles in this sector-time_bin.
    #Seperate into P.A bins of 10 degrees.                
    pa_sector = math.floor(pa*R2D/pa_bin)
    pa_flux[sector][timestep][pa_sector] += 1       #Number of detected particles in this sector-time_bin-pa_bin.


################################### CROSSING PARTICLES LAMDA-TIME PLOT #####################################
fig, ax = plt.subplots()

ax.scatter(detected_time, detected_lamda*R2D, c = detected_id)

ax.grid(alpha=0.8)

ax.set(xlabel="time(s)", ylabel="latitude(deg)", xlim=(0,1), xticks=np.arange(0,t+time_bin,time_bin))
plt.title("Detected particles in time. [Population: " +str(population)+ ", lamda0: " +str(round(lamda0[0]*R2D))+ "]\n"+"Northward particles are captured below the satellite.\nSouthward particles are captured above the satellite.")
plt.annotate("SATELLITE",xy=(0.45,30.0002),color="blue",weight="semibold")
ax.ticklabel_format(useOffset=False)    #disable e notation.
ax.axhline(y = 30 ,color="b", linestyle="dashed")

plt.savefig("simulation_MM/Crossing_particles.png", dpi=100)


######################################### PARTICLE SUM - 360 PLOT ###########################################
sum_flux = [0 for i in range(timesteps)]
time_min  = 0   #time_bin bounds
time_max  = time_bin
for timestep in range(0,timesteps):
    sum_flux[timestep] = sum(time_min<t_i<time_max for t_i in detected_time) #detected_time is sorted  
    time_min += time_bin
    time_max += time_bin


fig, ax = plt.subplots()

for timestep in range(0,timesteps):
    ax.set(xlabel="time", ylabel="sum", xticks=np.arange(0,t+time_bin,time_bin))
    ax.grid(alpha=0.8,axis="x")
    plt.title("Detected particle sum in all look_dirs")
    if sum_flux[timestep]!=0:
        ax.scatter(timestep*0.1+0.1/2,sum_flux[timestep], s=sum_flux[timestep])
plt.savefig("simulation_MM/Particle_sum.png", dpi=100)


#################################### PARTICLE COUNT - LOOK_DIRS MOVIE ########################################
fig, ax = plt.subplots(2,figsize=(7,7))

FFMpegWriter = manimation.writers["ffmpeg"]
metadata = dict(title="Particle bining",comment="Time bins of 0.1s, look directions of 15 degrees")
fps = 1
writer = FFMpegWriter(fps=fps, metadata = metadata)
print("Generating mp4 file...\nDuration of mp4 file will be:",(sectors*timesteps*fps)/60, "minutes")

with writer.saving(fig, "simulation_MM/Detected_particles.mp4", 100):

    for sector in range(0,sectors):

        for timestep in range(0,timesteps):           
            
            ax[0].set_xlabel(xlabel="time bins(s)")
            ax[0].set_ylabel(ylabel="count")
            ax[0].set_title("Particle count in Sector "+str(sector)+" - Look direction: ("+str(sector*sector_range-sector_range/2)+"\N{DEGREE SIGN},"+str(sector*sector_range+sector_range/2)+"\N{DEGREE SIGN})",loc="left",fontdict=font)              
            ax[1].set_xlabel("P.A bins(deg)") 
            ax[1].set_ylabel("count") 
            ax[1].set_title("P.A dstr in Timestep "+str(timestep), loc="left", fontdict=font,x=-0.015)
        
            for timestep2 in range(0,timesteps):           
                if sctr_flux[sector][timestep2] !=0: #to show only non zero dots.
                    ax[0].scatter(timestep2*time_bin+(time_bin/2),sctr_flux[sector][timestep2],c=colors[timestep2])

            for pa_sector in range(0,pa_sectors):     
                if pa_flux[sector][timestep][pa_sector]!=0: #to show only non zero dots.
                    ax[1].scatter(pa_sector*pa_bin+(pa_bin/2), pa_flux[sector][timestep][pa_sector],c=colors[timestep])  #y axis? flux? p.a can be very close which makes dots merge
            
            ax[0].grid(alpha=0.8, axis="x" )
            ax[1].grid(alpha=.8, axis="x") 
            ax[0].set_xticks(ticks=np.arange(0,t+time_bin,time_bin)) #time_bin ticks
            ax[0].set_xticks(ticks=np.arange(time_bin/2,t,time_bin),minor=True) #timestep ticks
            ax[0].set_xticklabels(labels=np.arange(0,timesteps,1),minor=True,c="red")  
            ax[1].set_xticks(ticks=np.arange(0,180+pa_bin,pa_bin)) #pa_bin ticks
            ax[1].set_xticklabels(labels=np.arange(0,180+pa_bin,pa_bin), rotation=90) #rotation 90
            ax[1].set_xticks(ticks=np.arange(pa_bin/2,180,pa_bin),minor=True) #pa_sector ticks
            ax[1].set_xticklabels(labels=np.arange(0,pa_sectors,1),minor=True, c="red") #pa_sector ticks
            #ax[1].set_yticks(ticks=np.linspace(0,max(pa_flux[sector][timestep]),10) ) 
            

            writer.grab_frame()
            ax[1].clear() #clear data of subplot 1 in every timestep.

        ax[0].clear() #clear data of subplot 2 in every sector change. 




























