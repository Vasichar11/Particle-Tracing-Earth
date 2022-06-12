import h5py
import numpy as np
import matplotlib.pyplot as plt
import csv 
import math
D2R=np.pi/180
R2D=1/D2R

view = 180 
sector_range = 1 #P.A bins #1deg
sectors = int(view/sector_range)

#Distribution read
f1 = h5py.File("h5files/10000p_normalAEQ_normalLAMDA_uniformETA_uniformEKIN.h5","r")
lamda0        = f1["lamda0"][()] #states when noWPI starts
ppar0         = f1["ppar0"][()]
pper0         = f1["pper0"][()]
alpha0        = f1["alpha0"][()]
aeq0          = f1["aeq0"][()]
eta0          = f1["eta0"][()]
Ekin0         = f1["Ekin0"][()]
time0         = f1["time0"][()]
aeq0_bins     = f1["aeq0_bins"][()]

f1.close()
#
#f2 = h5py.File("h5files/5000p_both.h5","r")
#population     = f2["population"][()]
#lamda00        = f2["lamda00"][()] #states when noWPI ends
#ppar00         = f2["ppar00"][()]
#pper00         = f2["pper00"][()]
#alpha00        = f2["alpha00"][()]
#aeq00          = f2["aeq00"][()]
#eta00          = f2["eta00"][()]
#time00         = f2["time00"][()]
#neg_id         = f2["neg_id"][()] #initial states of particles that develop negative p.a
#nan_id         = f2["nan_id"][()] #initial states of particles that develop nan
#high_id        = f2["high_id"][()] #initial states of particles that develop high p.a
#



############################################### SAVE CSV FILES ######################################################33
#Initials to CSV file
header = ['id', 'aeq0_deg', 'lamda0_deg', 'ppar0', 'pper0', 'alpha0_deg', 'eta0', 'Ekin(keV)']
data = []
id=0
for aq0,l0,par0,per0,a0,e0,E0 in zip(aeq0, lamda0, ppar0, pper0, alpha0, eta0, Ekin0):
    data.extend([[id,aq0*R2D,l0*R2D,par0,per0,a0*R2D,e0,Ekin0]])
    id=id+1 #id is the element's location in the list since the particles were saved like that
with open("dstr_data.csv", "w") as file1:
    writer = csv.writer(file1)
    writer.writerow(header)
    writer.writerows(data)



###################################### PLOT PIE INITIAL DISTRIBUTION ####################################

##BINNING
lamda_range = 10
aeq_range = 10
eta_range = 10
Ekin_range = 10


lamda_sectors = int((max(lamda0)-min(lamda0))*R2D / lamda_range)
aeq_sectors   = int((max(aeq0)-min(aeq0))*R2D / aeq_range)
eta_sectors   = int((max(eta0)-min(eta0))*R2D  / eta_range)
Ekin_sectors  = int((max(Ekin0)-min(Ekin0))*R2D / Ekin_range)
print(lamda_sectors, aeq_sectors, eta_sectors, Ekin_sectors)


def thread1_binning(lamda_range): 
    lamda_bins = [0 for i in range (lamda_sectors)]
    for lamda in lamda0: 
        lamda_bins[math.floor(lamda*R2D/lamda_range)] += 1  
    return lamda_bins


def thread2_binning(aeq_range): 
    aeq_bins   = [0 for i in range (aeq_sectors)]
    for aeq in aeq0: 
        aeq_bins  [math.floor(aeq/aeq_range)] += 1
    return aeq_bins


def thread3_binning(eta_range): 
    eta_bins   = [0 for i in range (eta_sectors)]
    for eta in eta0: 
        eta_bins  [math.floor(eta/eta_range)] += 1
    return eta_bins


def thread4_binning(Ekin_range): 
    Ekin_bins  = [0 for i in range (Ekin_sectors) ]
    for Ekin in Ekin0: 
        Ekin_bins [math.floor(Ekin/Ekin_range)] += 1
    return Ekin_bins

lamda_bins=thread1_binning(lamda_range)
aeq_bins=thread2_binning(aeq_range)
eta_bins=thread3_binning(eta_range)
Ekin_bins=thread4_binning(Ekin_range)


lamda_labels = np.arange(min(lamda0*R2D),max(lamda0*R2D),lamda_range) 
lamda_labels = [str(round(n,2)) for n in lamda_labels] 
aeq_labels   = np.arange(min(aeq0*R2D),max(aeq0*R2D),lamda_range)
aeq_labels   = [str(round(n,2)) for n in aeq_labels] 
eta_labels   = np.arange(min(eta0*R2D),max(eta0*R2D),lamda_range)
eta_labels   = [str(round(n,2)) for n in eta_labels] 
Ekin_labels  = np.arange(min(Ekin0),max(Ekin0),lamda_range)
Ekin_labels  = [str(round(n,2)) for n in Ekin_labels] 

# Make figure and axes
fig, axs = plt.subplots(2, 2)
print(lamda_labels)
print(lamda_bins)
print(aeq_labels)
print(eta_labels)
print(Ekin_labels)
# A standard pie plot
axs[0, 0].pie(lamda_bins, labels=lamda_labels, autopct='%1.1f%%', shadow=True)
axs[0, 1].pie(aeq_bins, labels=aeq_labels, autopct='%1.1f%%', shadow=True)
axs[1, 0].pie(eta_bins, labels=eta_labels, autopct='%1.1f%%', shadow=True)
axs[1, 1].pie(Ekin_bins, labels=Ekin_labels, autopct='%1.1f%%', shadow=True)

fig.savefig("simulation_MM/Distribution_Pie.png",dpi=200)



"""





######################################## PLOT INITIAL DISTRIBUTION #######################################

#Distribution data in one plot
fig, ax = plt.subplots(2,2)

#Ekin and eta distribution
ax[0,0].scatter(aeq0*R2D,Ekin0)
ax[0,0].grid(alpha=.3)
ax[0,0].set(xlabel="Aeq(deg)",ylabel="Ekin",title="Initial aeq-Ekin of simulated particles")
ax[1,0].scatter(aeq0*R2D,eta0)
ax[1,0].grid(alpha=.3)
ax[1,0].set(xlabel="Aeq(deg)",ylabel="eta",title="Initial aeq-eta of simulated particles")

#P.A distribution in bins
for sec in range(0,sectors):
    ax[0,1].scatter(sec,aeq0_bins[sec],s=2,alpha=1)
ax[0,1].grid(alpha=.3)
ax[0,1].set(xlabel="Sectors",xlim=(0,sectors-1),xticks=np.arange(0,sectors,10),ylabel="dN",title="Aeq0 distribution, sector range "+str(sector_range)+" degrees")
ax[0,1].set_yscale("log")


#P.A and latitude distribution
ax[1,1].scatter(lamda0*R2D,aeq0*R2D,s=0.5,alpha=0.1)
ax[1,1].grid(alpha=.3)
ax[1,1].set(xlabel="Latitude(deg)",ylabel="Equatorial P.A",title="Initial lat-aeq of simulated particles",ylim=(1,179),xlim=(-90,90),xticks=np.linspace(-90,90,5))
ax[1,1].axhline(y = 90, color ="b", linestyle="dashed")
fig.savefig("simulation_MM/Distribution_plot.png",dpi=200)

"""