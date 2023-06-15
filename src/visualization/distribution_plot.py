import h5py
import numpy as np
import matplotlib.pyplot as plt
import csv 
import math
import threading
D2R=np.pi/180
R2D=1/D2R

view = 180 
sector_range = 1 #P.A bins #1deg
sectors = int(view/sector_range)

#Distribution read
f1 = h5py.File("output/files/10000p_normalAEQ_normalLAMDA_constantETA_constantEKIN.h5","r")
lamda0        = f1["lamda0"][()] #states when noWPI starts
ppar0         = f1["ppar0"][()]
pper0         = f1["pper0"][()]
alpha0        = f1["alpha0"][()]
aeq0          = f1["aeq0"][()]
eta0          = f1["eta0"][()]
Ekin0         = f1["Ekin0"][()]
time0         = f1["time0"][()]
aeq0_bins     = f1["aeq0_bins"][()]
population    = len(lamda0) 
f1.close()
#
#f2 = h5py.File("output/files/5000p_both.h5","r")
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
    data.extend([[id,aq0*R2D,l0*R2D,par0,per0,a0*R2D,e0,E0]])
    id=id+1 #id is the element's location in the list since the particles were saved like that
with open("dstr_data.csv", "w") as file1:
    writer = csv.writer(file1)
    writer.writerow(header)
    writer.writerows(data)



###################################### PLOT PIE INITIAL DISTRIBUTION ####################################

##BINNING
lamda_range   = 5 #Bins of 10 degrees
lamda_domain  = (max(lamda0)-min(lamda0))*R2D
lamda_sectors = math.ceil(lamda_domain/lamda_range) #Ceil because we have also the "0" bin.
lamda_bins    = [0 for i in range (lamda_sectors)]
lamda_labels  = np.arange(min(lamda0*R2D),max(lamda0*R2D)+  lamda_range,lamda_range) #np.arange doesn't include endpoint
lamda_labels2 = []
for i in range(len(lamda_labels)-1):
    temp_tuple = (str(int(lamda_labels[i])),str(int(lamda_labels[i+1])) )
    lamda_labels2.append(" to ".join(temp_tuple)) #Round labels - make them ranges

aeq_range    = 20
aeq_domain   = (max(aeq0)-min(aeq0))*R2D
aeq_sectors  = math.ceil(aeq_domain / aeq_range)
aeq_bins     = [0 for i in range (aeq_sectors)]
aeq_labels   = np.arange(min(aeq0*R2D),max(aeq0*R2D),aeq_range)
aeq_labels   = [str(int(n)) for n in aeq_labels] 
aeq_labels2  = []
for i in range(len(aeq_labels)-1):
    temp_tuple = (str(int(aeq_labels[i])),str(int(aeq_labels[i+1])) )
    aeq_labels2.append(" to ".join(temp_tuple)) 

eta_range    = 40
eta_domain   = (max(eta0)-min(eta0))*R2D
eta_sectors  = math.ceil(eta_domain / eta_range)
eta_bins     = [0 for i in range (eta_sectors)]
eta_labels   = np.arange(min(eta0*R2D),max(eta0*R2D),eta_range)
eta_labels   = [str(int(n)) for n in eta_labels] 
eta_labels2  = []
for i in range(len(eta_labels)-1):
    temp_tuple = (str(int(eta_labels[i])),str(int(eta_labels[i+1])) )
    eta_labels2.append(" to ".join(temp_tuple)) 

Ekin_range   = 100 #Bins of 100kev 
Ekin_domain  = (max(Ekin0)-min(Ekin0))
Ekin_sectors = math.ceil(Ekin_domain / Ekin_range)
Ekin_bins    = [0 for i in range (Ekin_sectors)]
Ekin_labels  = np.arange(min(Ekin0),max(Ekin0),Ekin_range)
Ekin_labels  = [str(int(n)) for n in Ekin_labels] 
Ekin_labels2 = []
for i in range(len(Ekin_labels)-1):
    temp_tuple = (str(int(Ekin_labels[i])),str(int(Ekin_labels[i+1])) )
    Ekin_labels2.append(" to ".join(temp_tuple)) 




def thread1_binning(): 
    for lamda in lamda0: 
        lamda_bins[math.floor((lamda-min(lamda0))*R2D/lamda_range)] += 1  
def thread2_binning(): 
    for aeq in aeq0: 
        aeq_bins[math.floor((aeq-min(aeq0))*R2D/aeq_range)] += 1  
def thread3_binning(): 
    for eta in eta0: 
        eta_bins  [math.floor((eta-min(eta0))*R2D/eta_range)] += 1
def thread4_binning(): 
    for Ekin in Ekin0: 
        Ekin_bins  [math.floor((Ekin-min(Ekin0))/Ekin_range)] += 1




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


fig, axs = plt.subplots(2, 2)
fig.suptitle('Particle Distribution Pies',weight="bold")

axs[0,0].pie(lamda_bins)
fig.legend(labels=lamda_labels2,loc="upper left")
axs[0,0].set(title="lamda dstr(deg)")

axs[0,1].pie(aeq_bins)
fig.legend(labels=aeq_labels2,loc="upper right")
axs[0,1].set(title="aeq dstr(deg)")

axs[1,0].pie(eta_bins)
fig.legend(labels=eta_labels2,loc="lower left")
axs[1,0].set(title="eta dstr(deg)")

axs[1,1].pie(Ekin_bins)
fig.legend(labels=Ekin_labels2,loc="lower right", prop={'size': 8})
axs[1,1].set(title="Ekin dstr(keV)")

fig.savefig("output/plots/"+str(population)+"p_Distribution_Pies.png",dpi=200)



######################################## PLOT INITIAL DISTRIBUTION #######################################
fig, ax = plt.subplots(2)

#P.A distribution in bins
for sec in range(0,sectors):
    ax[0].scatter(sec,aeq0_bins[sec],s=2,alpha=1)
ax[0].grid(alpha=.3)
ax[0].set(xlim=(0,sectors-1),xticks=np.arange(0,sectors,10),ylabel="dN",title="Sector range "+str(sector_range)+" deg")
ax[0].set_yscale("log")

#P.A and latitude distribution
ax[1].scatter(lamda0*R2D,aeq0*R2D,s=0.5,alpha=0.1)
ax[1].grid(alpha=.3)
ax[1].set(xlabel="Latitude(deg)",ylabel="Equatorial P.A",ylim=(1,179),xlim=(-90,90),xticks=np.linspace(-90,90,5))
ax[1].axhline(y = 90, color ="b", linestyle="dashed")

fig.suptitle('Particle aeq-lamda distributions',weight="bold")
fig.savefig("output/plots/"+str(population)+"p_Distribution_Plots.png",dpi=200)
