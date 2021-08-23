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
f2.close();
detected_pop   = len(detected_lamda)
##################################### TELESCOPE SPECIFICATION #########################################
t = round(time_sim[0][-1]) #simulation time
tRes = 0.1 #seconds to distinquish events
timesteps = int (t / tRes) 
view = 360
sector_range = 15
sectors = int(view/sector_range)
########################################### FONT DICT  ################################################
font = {'family': 'serif',
        'color':  'blue',
        'weight': 'medium',
        'size': 10,
        }
################################# CROSSING PARTICLE LAMDA-TIME PLOT ####################################
fig, ax = plt.subplots()

ax.scatter(detected_time, detected_lamda*R2D, c = detected_id)

ax.grid(alpha=0.8)

ax.set(xlabel="time(s)", ylabel="latitude(deg)", xlim=(0,1), xticks=np.arange(0,t+tRes,tRes))
plt.title("Detected particles in time. [Population: " +str(population)+ ", lamda0: " +str(round(lamda0[0]*R2D))+ "]\n"+"North heading particles are captured below the satellite.\nSouth heading particles are captured above the satellite.")
plt.annotate("SATELLITE",xy=(0.45,30.0002),color="blue",weight="semibold")
ax.ticklabel_format(useOffset=False)    #disable e notation.
ax.axhline(y = 30 ,color="b", linestyle="dashed")

plt.show()
################################ Sorting: Particles in Temporal order ####################################
r=0
while r < len(detected_time) - 1:
    
    if(detected_time[r+1]<detected_time[r]):
        
        detected_time[r], detected_time[r+1]  = detected_time[r+1], detected_time[r]   #swap
        detected_lamda[r],detected_lamda[r+1] = detected_lamda[r+1],detected_lamda[r]  
        detected_id[r],   detected_id[r+1]    = detected_id[r+1],   detected_id[r]  
        
        r = -1 #next loop r=0 => sort from start.
    
    r = r+1
#print(detected_time)





##########################################################################################################
######################################### POST PROCESSING ################################################
##########################################################################################################





####################################### PARTICLE SUM - 360 PLOT #########################################
flux = [0 for i in range(timesteps)]
time_min  = 0   #time_bin bounds
time_max  = tRes
for timestep in range(0,timesteps):
    flux[timestep] = sum(time_min<i<time_max for i in detected_time) #detected_time is sorted
    time_min += tRes
    time_max += tRes
#print(flux)

fig, ax = plt.subplots()

for timestep in range(0,timesteps):
    ax.set(xlabel="time", ylabel="sum", xticks=np.arange(0,t+tRes,tRes))
    ax.grid(alpha=0.8,axis="x")
    plt.title("Detected particle sum in all look_dirs")
    if flux[timestep]!=0:
        ax.scatter(timestep*0.1+0.1/2,flux[timestep], s=flux[timestep])
plt.show()




############################### ARRANGE PARTICLE FLUX-P.A DSTR ARRAYS ######################################
flux2=[ [0 for i in range(timesteps)] for j in range(sectors) ]   #flux2[timesteps][sectors]
time_min  = 0   #reset
time_max  = tRes
pa_dstr=[]
pa_dstr_time=[]


for i in range(0,sectors):      #2d lists, with a nested list to store P.A distribution of detected particles and their detection time.
    sect = []
    for j in range(0,timesteps):   
        sect.append([])         #Append 10 empty lists in each sector. In these empty lists, P.A dstr will be placed.
    pa_dstr.append(sect)

for i in range(0,sectors):     
    sect2 = []
    for j in range(0,timesteps):   
        sect2.append([])        
    pa_dstr_time.append(sect2)

for timestep in range(0,timesteps,1):    
    sector=0
    for look_dir in range(0,view,sector_range):
        for time,pa in zip(detected_time,detected_alpha):   #iterate in both arrays.
            if ( time_min<time<time_max  and  (int(pa*R2D+180)%360) >= (look_dir-sector_range/2)  and  (int(pa*R2D+180)%360) <= (look_dir+sector_range/2)) :
                #print("Particle with",time,pa*R2D,"is detected when",timestep,look_dir)
                flux2[sector][timestep] += 1             #Number of detected particles in this sctor-timestep.
                pa_dstr[sector][timestep].append(pa*R2D) #Their p.a dstr.
                pa_dstr_time[sector][timestep].append(time) #time when these are detected.
        sector += 1     #Next sector/look_dir
    time_min += tRes    #Next time_bin
    time_max += tRes
#Particles are sorted in temporal order, it's not efficient to iterate through all particles in the nested loop.
#1: don't sort or 2: use while loop.
#print(flux2)

#for i in range(0,detected_pop):
#    print(detected_time[i],"s",detected_alpha[i]*R2D,"\N{DEGREE SIGN}","detection look_dir: ",int(detected_alpha[i]*R2D + 180)%360)






########################################### MOVIE DATA ###################################################
FFMpegWriter = manimation.writers["ffmpeg"]
metadata = dict(title="Particle bining",comment="Time bins of 0.1s, look directions of 15 degrees")
fps = 1
writer = FFMpegWriter(fps=fps, metadata = metadata)




################################# PARTICLE COUNT - LOOK_DIRS PLOT ########################################
print("Generating mp4 file...\nDuration of mp4 file will be:",(sectors*timesteps*fps)/60, "minutes")

fig, ax = plt.subplots(2,figsize=(7,7))
colors = np.array(["red","green","blue","yellow","pink","black","orange","purple","beige","brown"]) #to seperate counts in timestep

with writer.saving(fig, "detected_particles.mp4", 100):
    
    for sector in range(13,16):         # sectors x timesteps = 24 x 10 = 240 plots --> movie?
    
        for timestep in range(0,timesteps):
    
            for timestep2 in range(0,timesteps):
    
                ax[0].set_xlabel(xlabel="time(s) in black and timestep in red")
                ax[0].set_ylabel(ylabel="count")
                ax[0].grid(alpha=0.8, axis="x" )
                ax[0].set_title("Particle count in Sector "+str(sector)+": ("+str(sector*sector_range-sector_range/2)+"\N{DEGREE SIGN},"+str(sector*sector_range+sector_range/2)+"\N{DEGREE SIGN})",loc="left",fontdict=font)    
                ax[0].set_xticks(ticks=np.arange(0,t+tRes,tRes)) #timestep ticks
                ax[0].set_xticks(ticks=np.arange(tRes/2,t,tRes),minor=True) #timestep ticks
                ax[0].set_xticklabels(labels=np.arange(0,10,1),minor=True,c="red") 
    
                if flux2[sector][timestep2] !=0: #to show only non zero dots.
                    ax[0].scatter(timestep2*0.1+(0.1/2),flux2[sector][timestep2],c=colors[timestep2])    
            
            #after finishing scattering for the first subplot
    
            ax[1].set_xlabel("detection time(s)") 
            ax[1].set_ylabel("P.A(degrees)") 
            ax[1].set_title("P.A dstr in Timestep "+str(timestep), loc="left", fontdict=font, y=-0.2)
            ax[1].grid(alpha=.8, axis="x")  
            ax[1].ticklabel_format(useOffset=False)  
    
            for p in range(0,flux2[sector][timestep]):  
                ax[1].scatter(pa_dstr_time[sector][timestep][p], pa_dstr[sector][timestep][p],c=colors[timestep])  #y axis? flux? p.a can be very close which makes dots merge
    
            writer.grab_frame()
            ax[1].clear() #clear data of subplot 1 in every timestep.

        ax[0].clear() #clear data of subplot 2 in every sector change. 

