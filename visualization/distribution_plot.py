import h5py
import numpy as np
import matplotlib.pyplot as plt
import csv 
import math
import threading
import os
D2R=np.pi/180
R2D=1/D2R

# Parameters
aeq_sector_deg = 1
latitude_sector_deg = 10
eta_sector_deg = 40
Ekin_sector = 100

# Read particle distribution from file entered by user
files_dir = os.path.join("output","files")  # Path to the output directory
h5_files = [file for file in os.listdir(files_dir) if (file.startswith("dstr") and file.endswith('.h5'))]
plots_dir = os.path.join("output","plots")


# Display the list of .h5 files
print("Available particle distribution .h5 files:")
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
    else:
        print("Invalid selection.")
except ValueError:
    print("Invalid input. Please enter a number.")

distribution_name = os.path.basename(selected_file).split(".h5")[0]

# Read file
f1 = h5py.File(selected_file,"r")
latitude0 = f1["latitude0"][()] 
ppar0 = f1["ppar0"][()]
pper0 = f1["pper0"][()]
alpha0 = f1["alpha0"][()]
aeq0 = f1["aeq0"][()]
eta0 = f1["eta0"][()]
Ekin0 = f1["Ekin0"][()]
time0 = f1["time0"][()]
population = len(latitude0) 
f1.close()
latitude0_deg = latitude0 * R2D
aeq0_deg = aeq0 * R2D
eta0_deg = eta0 * R2D

# Sectors
aeq_domain  = 180
latitude_domain  = 180
eta_domain  = 360
Ekin_domain  = np.ptp(Ekin0) 
aeq_sectors = math.ceil(aeq_domain/aeq_sector_deg)
latitude_sectors = math.ceil(latitude_domain/latitude_sector_deg)
eta_sectors = math.ceil(eta_domain/eta_sector_deg)
Ekin_sectors = math.ceil(Ekin_domain / Ekin_sector)

############################################### SAVE CSV FILE ###########################################

# Initials to CSV file
header = ['id', 'aeq0_deg', 'latitude0_deg', 'ppar0', 'pper0', 'alpha0_deg', 'eta0', 'Ekin(keV)']
data = []
id=0
for aq0,l0,par0,per0,a0,e0,E0 in zip(aeq0, latitude0, ppar0, pper0, alpha0, eta0, Ekin0):
    data.extend([[id,aq0*R2D,l0*R2D,par0,per0,a0*R2D,e0*R2D,E0]])
    id=id+1 #id is the element's location in the list since the particles were saved like that
h5_file = selected_file.split()[0].split(".h5")
csv_file = os.path.join(files_dir, h5_files[selection-1].split("h5")[0] +"csv")
with open(csv_file, "w") as file1:
    writer = csv.writer(file1)
    writer.writerow(header)
    writer.writerows(data)

###################################### BINNING BASED ON SECTOR RANGES ####################################

latitude_bins = [0] * latitude_sectors
latitude_labels = np.linspace(min(latitude0_deg), max(latitude0_deg), latitude_sectors)
latitude_labels2 = [f"{int(start)} to {int(end)}" for start, end in zip(latitude_labels[:-1], latitude_labels[1:])]

aeq_bins = [0] * aeq_sectors
aeq_labels = np.linspace(min(aeq0_deg), max(aeq0_deg), aeq_sectors)
aeq_labels2 = [f"{int(start)} to {int(end)}" for start, end in zip(aeq_labels[:-1], aeq_labels[1:])]

eta_bins = [0] * eta_sectors
eta_labels = np.linspace(min(eta0_deg), max(eta0_deg), eta_sectors)
eta_labels2 = [f"{int(start)} to {int(end)}" for start, end in zip(eta_labels[:-1], eta_labels[1:])]

Ekin_bins = [0] * Ekin_sectors
Ekin_labels = np.linspace(min(Ekin0), max(Ekin0), Ekin_sectors)
Ekin_labels2 = [f"{int(start)} to {int(end)}" for start, end in zip(Ekin_labels[:-1], Ekin_labels[1:])]

###################################### ONE THREAD PER BINNING PROCESS ####################################

def thread1_binning(): 
    for latitude_deg in latitude0_deg:
        latitude_bins[math.floor(latitude_deg/latitude_sector_deg)-1] += 1  
def thread2_binning(): 
    for aeq_deg in aeq0_deg: 
        aeq_bins[math.floor(aeq_deg/aeq_sector_deg)-1] += 1  
def thread3_binning(): 
    for eta_deg in eta0_deg: 
        eta_bins[math.floor(eta_deg/eta_sector_deg)-1] += 1
def thread4_binning(): 
    for Ekin in Ekin0: 
        Ekin_bins[math.floor(Ekin/Ekin_sector)-1] += 1

th1 = threading.Thread(target=thread1_binning)
th2 = threading.Thread(target=thread2_binning)
th3 = threading.Thread(target=thread3_binning)
th4 = threading.Thread(target=thread4_binning)

th1.start()
th2.start()
th3.start()
th4.start()
th1.join()
th2.join()
th3.join()
th4.join()

###################################### PLOT PIES OF INITIAL DISTRIBUTION ####################################

fig, axs = plt.subplots(2, 2)
fig.suptitle('Particle Distribution Pies',weight="bold")
axs[0,0].pie(latitude_bins)
fig.legend(labels=latitude_labels2,loc="upper left")
axs[0,0].set(title="latitude dstr(deg)")

axs[0,1].pie(aeq_bins)
fig.legend(labels=aeq_labels2,loc="upper right")
axs[0,1].set(title="aeq dstr(deg)")

axs[1,0].pie(eta_bins)
fig.legend(labels=eta_labels2,loc="lower left")
axs[1,0].set(title="eta dstr(deg)")

axs[1,1].pie(Ekin_bins)
fig.legend(labels=Ekin_labels2,loc="lower right", prop={'size': 8})
axs[1,1].set(title="Ekin dstr(keV)")

fig.savefig(os.path.join(plots_dir, distribution_name + "p_pies.png"),dpi=200)

######################################## PLOT INITIAL DISTRIBUTION #######################################

fig, ax = plt.subplots(2)
for sec in range(0,aeq_sectors):
    if aeq_bins[sec] > 0:
        ax[0].scatter(sec,aeq_bins[sec],s=2,alpha=1, color="blue")
ax[0].grid(alpha=0.3)
ax[0].set(ylabel="dN", title="Sector range " + str(aeq_sector_deg) + " deg", xlim=(np.min(aeq0_deg), np.max(aeq0_deg)))
ax[0].axvline(x = 90, color ="r", linestyle="dashed")
ax[0].set_yscale("log")

# P.A and latitude distribution
ax[1].scatter(latitude0*R2D,aeq0*R2D,s=0.5,alpha=0.1, color="blue")
ax[1].grid(alpha=.3)
ax[1].set(xlabel="Latitude(deg)",ylabel="Equatorial P.A",ylim=(np.min(aeq0_deg),np.max(aeq0_deg)),xlim=(-90,90),xticks=np.linspace(-90,90,5))
ax[1].axhline(y = 90, color ="r", linestyle="dashed")

fig.suptitle('Particle aeq-latitude distributions',weight="bold")
fig.savefig(os.path.join(plots_dir, distribution_name + "p_aeq_lamda.png"),dpi=200)

print("Plots saved at: ",plots_dir)