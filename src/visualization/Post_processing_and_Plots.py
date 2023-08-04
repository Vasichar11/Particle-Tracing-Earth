import numpy as np 
import sys
import matplotlib.pyplot as plt
import h5py
import math
import matplotlib.animation as manimation
from random import randint
import threading
import os

############################################### CONSTANTS ###############################################

np.set_printoptions(threshold=sys.maxsize)
D2R=np.pi/180
R2D=1/D2R
files_dir = os.path.join("output","files")  # Filepath for files
plots_dir = os.path.join("output","plots")  # Filepath for plots

############################################# READ HDF5 ###################################################

# Function to read particle distribution from file entered by user
def get_file_selection(h5_files, files_dir):
    # Display the list of .h5 files
    print("Available simulation .h5 files:")
    for i, file in enumerate(h5_files):
        print(f"{i+1}. {file}")
    
    # Prompt the user to select a file
    selection = input("Enter the number corresponding to the file you want to select: ")
    
    # Validate the user's input and retrieve the selected filepath
    try:
        selection = int(selection)
        if 1 <= selection <= len(h5_files):
            selected_file = os.path.join(files_dir, h5_files[selection-1])
            print(f"You selected: {selected_file}")
            return selected_file
        else:
            print("Invalid selection.")
    except ValueError:
        print("Invalid input. Please enter a number.")

# Get the list of Simulation1 .h5 files
print("Select Simulation1 simulation")
h5_files = [file for file in os.listdir(files_dir) if (file.startswith("sim") and file.endswith('.h5'))]

# Prompt the user for the first simulation file selection
first_file = get_file_selection(h5_files, files_dir)

# Read file
f1 = h5py.File(first_file,"r")
detected_latitude1 = f1["ODPT.latitude"][()]
detected_time1 = f1["ODPT.time"][()]
detected_id1 = f1["ODPT.id"][()]
detected_alpha1 = f1["ODPT.alpha"][()]
detected_aeq1 = f1["ODPT.aeq"][()]
telescope_latitude1 = f1["ODPT.latitude_deg"][()]
population1 = f1["population"][()]
t1 = f1["t"][()]
precip_id1 = f1["precip_id"][()]
precip_latitude1 = f1["precip_latitude"][()]
precip_aeq1 = f1["precip_aeq"][()]
precip_time1 = f1["precip_time"][()]
neg_id1 = f1["neg_id"][()]
high_id1 = f1["high_id"][()]
nan_id1 = f1["nan_id"][()]
f1.close()

# Particles that should be excluded 
excluded_ids = []
excluded_ids.extend(neg_id1)
excluded_ids.extend(nan_id1)
excluded_ids.extend(high_id1)

# Prompt the user for the second simulation file selection (optional)
read_second_file = None
comparing_simulations = False
while read_second_file!="yes" and read_second_file!="no":

    read_second_file = input("Do you want to compare this simulation with another simulation? (yes/no): ")
    
    if read_second_file.lower() == "yes":
        second_file = get_file_selection(h5_files, files_dir)
        comparing_simulations = True
        break
    elif read_second_file.lower() == "no":
        comparing_simulations = False
        break
    else:
        print("Invalid input. Please enter 'yes' or 'no'.")


if comparing_simulations:
    f2 = h5py.File(second_file,"r")
    t2 = f2["t"][()]
    detected_latitude2 = f2["ODPT.latitude"][()]
    detected_time2 = f2["ODPT.time"][()]
    detected_id2 = f2["ODPT.id"][()]
    detected_alpha2 = f2["ODPT.alpha"][()]
    detected_aeq2 = f2["ODPT.aeq"][()]
    telescope_latitude2= f2["ODPT.latitude_deg"][()]
    population2 = f2["population"][()]
    precip_id2 = f2["precip_id"][()]
    precip_latitude2 = f2["precip_latitude"][()]
    precip_aeq2 = f2["precip_aeq"][()]
    precip_time2 = f2["precip_time"][()]
    neg_id2 = f2["neg_id"][()]
    high_id2 = f2["high_id"][()]
    nan_id2 = f2["nan_id"][()]
    f2.close()

    # Ensure that selected simulations are for same configuration
    if t1!=t2:
        raise Exception("Selected simulation files are not executed for the same duration. Simulation 1: %ds, Simulation 2: %ds" %(t1, t2))
    if telescope_latitude1!=telescope_latitude2:
        raise Exception("Selected simulation files are not executed with the satellite in the same location. Simulation 1: %ds, Simulation 2: %ds" %(telescope_latitude1, telescope_latitude2))
    if population1!=population2:
        raise Exception("Selected simulation files are not executed for the particle population. Simulation 1: %ds, Simulation 2: %ds" %(population1, population2))

    # Exclude the following particles as well
    excluded_ids.extend(neg_id2)
    excluded_ids.extend(nan_id2)
    excluded_ids.extend(high_id2)


############################ TELESCOPE SPECIFICATION -- BINNING PARAMETERS ################################

time_bin  = 0.1                # seconds to distinquish events(time resolution)
timesteps = math.ceil(t1 / time_bin) # t stops at the last timestep (e.g 14.9)
view = 180
sector_range = 2 # P.A bins, look directions
sectors = int(view/sector_range)



###################################### POST PROCESSING - PLOTS ###########################################

print()
print("Population     :",population1)
print("Simulation time:",t1,"seconds")
print("Time bins      :",time_bin,"   seconds") # Sampling frequency should comply with Nyquist
print("Timesteps      :",timesteps) 
print("P.A bins       :",sector_range,"degrees")
print("Sectors        :",sectors)
print()
print("|Sim1   Simu2|")
print(" ",len(precip_id1),"     ",len(precip_id2)," particles escaped")
print()
print("|Sim1   Sim2|")
print(" ",len(neg_id1),"     ",len(neg_id2),"  particles developed negative P.A")
print(" ",len(high_id1),"     ",len(high_id2),"  particles developed high P.A")
print(" ",len(nan_id1),"     ",len(nan_id2),"  particles developed nan P.A") 
print("These particles will be excluded from the bining population losing:", ((len(neg_id1)+len(neg_id2)+len(nan_id1)+len(nan_id2)+len(high_id1)+len(high_id2))/population1)*100,"% of the population\n")

########################################### FONTS AND COLORS #############################################
font = {'family': 'serif',
        'color':  'blue',
        'weight': 'bold',
        'size': 2}
colors = []
for i in range(max(timesteps,sectors)):      # Colors to seperate timesteps or sectors.
    colors.append('#%06X' % randint(0, 0xFFFFFF))
################################### CROSSING PARTICLES LAMDA-TIME PLOT ####################################
#Simulation1
#"""

# Find fixed limits for y axis:
y_max = max(np.max(detected_latitude1), np.max(detected_latitude2)) * R2D
y_min = min(np.min(detected_latitude1), np.min(detected_latitude2)) * R2D

fig, ax = plt.subplots()
ax.scatter(detected_time1, detected_latitude1*R2D, c = detected_id1, s=0.3, facecolors='none', edgecolors='r')
ax.grid(alpha=0.8)
ax.set(xlabel="time(s)", ylabel="latitude(deg)", ylim=(y_min,y_max))
plt.title("$Population$: " +str(population1)+"\nNorthward particles are captured below the satellite.\nSouthward particles are captured above the satellite",size="medium")
plt.annotate("SATELLITE",xy=(t1/2,telescope_latitude1+0.0002),color="blue",weight="semibold")
ax.ticklabel_format(useOffset=False)  # Disable e notation.
ax.axhline(y = telescope_latitude1 ,color="b", linestyle="dashed")
plt.savefig(plots_dir+"/Crossing_particles.png", dpi=100)
if comparing_simulations:
    #BOTH
    fig, ax = plt.subplots()
    ax.scatter(detected_time2, detected_latitude2*R2D, c = detected_id2, s=0.3, facecolors='none', edgecolors='r')
    ax.grid(alpha=0.8)
    ax.set(xlabel="time(s)", ylabel="latitude(deg)", ylim=(y_min,y_max))
    plt.title("$Population$: " +str(population2)+"\nNorthward particles are captured below the satellite.\nSouthward particles are captured above the satellite",size="medium")
    plt.annotate("SATELLITE",xy=(t2/2,telescope_latitude2+0.0002),color="blue",weight="semibold")
    ax.ticklabel_format(useOffset=False)  # Disable e notation.
    ax.axhline(y = telescope_latitude2 , color="b", linestyle="dashed")
    plt.savefig(plots_dir+"/Crossing_particles2.png", dpi=100)
#"""

"""
########################################## MINIMUM DETA_DT ################################################
#For every particle in the population the state for minimum deta_dt is saved. 
fix,ax = plt.subplots()
cc = ax.scatter(np.rad2deg(saved_latitude),saved_deta_dt,s=1, c=saved_Ekin,cmap='jet')
ax.axhline(y = 0 ,color="b", linestyle="dashed")
ax.set(xlabel="latitude[deg]",ylabel="deta_dt")
cbar=plt.colorbar(cc)
cbar.ax.set_title('Ekin[keV]')
plt.savefig(plots_dir+"/min_deta_dt.png")
"""

############################################## BINNING ####################################################
###BINNING WITH LOCAL PITCH ANGLE: detected_alpha1 not detected_aeq1
# Formula for aeq is not valid. Do the binning with local P.A since satellite is @0deg aeq~=alpha.
#"""
sctr_flux = np.zeros((sectors, timesteps)) #sctr_flux[sectors][timesteps]
sum_flux = np.zeros(timesteps)
sctr_flux2 = np.zeros((sectors, timesteps))
sum_flux2 = np.zeros(timesteps)
sctr_flux_precip = np.zeros((sectors, timesteps))
sum_flux_precip = np.zeros(timesteps)

def thread1_binning():
# Simulation1
    for time,pa,id in zip(detected_time1,detected_alpha1,detected_id1): # Iterate in both array elements. No sorting required.
        timestep = math.floor(time/time_bin)
        sector   = math.floor(pa*R2D/sector_range)
        if(sector>=sectors or timestep>=timesteps): # Uneeded? Happens when NAN particles are kept in the simulation(jumping steps) or when initializing close to 0 and 180 aeq. WHY?
            raise Exception("Simulation1 Particle",id,"at time",time,"needs sector",sector,"in timestep",timestep)
        if id in excluded_ids:
            continue # don't include these particles in the binning process
        sctr_flux[sector][timestep] += 1              # Number of detected particles in this sector-timestep.
        sum_flux[timestep] += 1                       # Sum of detected particles in this timestep from all look directions.
    print("Binning Simulation1 done!")

def thread2_binning():
# Simulation1 and then Simulation2
    for time,pa,id in zip(detected_time2,detected_alpha2,detected_id2):
        timestep = math.floor(time/time_bin)
        sector   = math.floor(pa*R2D/sector_range)
        if(sector>=sectors or timestep>=timesteps): 
            raise Exception("Simulation2 Particle",id,"with pa",pa*R2D,"at time",time,"needs sector",sector,"in timestep",timestep)
        # Do not include these particles
        if id in excluded_ids:
            continue
        sctr_flux2[sector][timestep] += 1              
        sum_flux2[timestep] += 1
    print("Binning Simulation1+Simulation2 done!")

def thread_precip_binning():
# Prepitating particles. Just iterating for the latter simulation(Simulation1 and then Simulation2) is enough.
# Particles cannot precipitate in Simulation1 simulation. They can only get lost before the first bounce.
# These particles escape with alpha close to 90(?), pitch angle binning to compare them with the particles that are crossing the equator
    for time,pa in zip(precip_time2,precip_aeq2):
        timestep = math.floor(time/time_bin)
        sector   = math.floor(pa*R2D/sector_range)
        if(sector>=sectors or timestep>=timesteps): 
            raise Exception("Precipitating Particle at time",time,"needs sector",sector,"in timestep",timestep)
        if id in excluded_ids:
            continue
        sctr_flux_precip[sector][timestep] += 1           
        sum_flux_precip[timestep] += 1
    print("Binning Precipitating done!")

#Create thread classes
th1=threading.Thread(target=thread1_binning)
th2=threading.Thread(target=thread2_binning)
th3=threading.Thread(target=thread_precip_binning)

th1.start()
if comparing_simulations:
    th2.start()
    th3.start()

th1.join()
th2.join()
th3.join()
# main waits until all threads are done

######################################### PARTICLE SUM - 360 PLOT ###########################################
#"""
fig, ax = plt.subplots()
plt.title("Detected particle sum in all look_dirs for "+str(t1)+" seconds, in "+str(timesteps)+" timesteps\n Satellite @"+str(telescope_latitude1)+" deg")
ax.set(xlabel="Time(s), in time_bins of "+str(time_bin)+"(s)", ylabel="Total Flux")
#ax.set_xticks(ticks=np.arange(0,t+time_bin,time_bin),minor=True) #ticks for time_bin seperation
ax.xaxis.grid(True, which='both')
for timestep in range(0,timesteps):
    if sum_flux[timestep]!=0:
        ax.scatter(timestep*time_bin+time_bin/2,sum_flux[timestep],s=10, c="black",label="Simulation1")
    if sum_flux2[timestep]!=0:  
        ax.scatter(timestep*time_bin+time_bin/2,sum_flux2[timestep],s=10, c="red", label="Simulation1 && Simulation2")  
    if timestep==0:# plot legend once
        ax.legend()
plt.savefig(plots_dir+"/Particle_sum_bins"+str(time_bin)+"s.png", dpi=100)
#"""
##################################### Simulation2-Simulation1 DIFF FOR HISTOGRAM #########################################
#"""
moved = np.zeros((timesteps, sectors))
precip = np.zeros((timesteps, sectors))
lost = np.zeros(sectors)

for sector in range(0,sectors):  
    for timestep in range(0,timesteps): 
        precip[timestep][sector] = sctr_flux_precip[sector][timestep] # Particles lost in this timestep and sector.
        lost[sector] += precip[timestep][sector]                      # Particles lost until this timestep,for this sector.
        moved[timestep][sector]  = abs(sctr_flux[sector][timestep] - sctr_flux2[sector][timestep]) - lost[sector] # Find difference between simulations and remove lost particles until this timestep for this sector
#"""       
###################################### (FLUX-P.A)*TIMESTEPS MOVIE ##########################################
#"""
fig,ax = plt.subplots()
FFMpegWriter = manimation.writers["ffmpeg"]
metadata2 = dict(title="P.A binning",comment="P.A bins of"+str(sector_range)+"degrees")
fps = 1
writer = FFMpegWriter(fps=fps, metadata = metadata2)
print("Generating P.A binning mp4 file...\nDuration of mp4 file will be:",(timesteps*fps), "seconds")
with writer.saving(fig, plots_dir+"/Bins_"+str(sector_range)+"deg_"+str(time_bin)+"s.mp4", 100):
    for timestep in range(0,timesteps):           
            
        for sector in range(0,sectors): 

            if(sctr_flux[sector][timestep]!=0): #plot only non zero dots
                ax.scatter(sector+0.5, sctr_flux[sector][timestep],c="black",s=1) 
            if(sctr_flux2[sector][timestep]!=0 and sctr_flux2[sector][timestep]!=sctr_flux[sector][timestep]): #plot dot only when it differs with Simulation1
                ax.scatter(sector+0.5, sctr_flux2[sector][timestep],c="red",s=1) 
             
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
