import h5py
import numpy as np
import matplotlib.pyplot as plt
import csv 

D2R=np.pi/180
R2D=1/D2R

view = 180 
sector_range = 10 #P.A bins #1deg
sectors = int(view/sector_range)

#Distribution read
f1 = h5py.File("h5files/10p_normalAEQ_normalLAMDA.h5","r")
lamda0        = f1["lamda0"][()] #states when noWPI starts
ppar0         = f1["ppar0"][()]
pper0         = f1["pper0"][()]
alpha0        = f1["alpha0"][()]
aeq0          = f1["aeq0"][()]
eta0          = f1["eta0"][()]
time0         = f1["time0"][()]
aeq0_bins     = f1["aeq0_bins"][()]
f1.close()

f2 = h5py.File("h5files/10p_both.h5","r")
population     = f2["population"][()]
lamda00        = f2["lamda00"][()] #states when noWPI ends
ppar00         = f2["ppar00"][()]
pper00         = f2["pper00"][()]
alpha00        = f2["alpha00"][()]
aeq00          = f2["aeq00"][()]
eta00          = f2["eta00"][()]
time00         = f2["time00"][()]
neg_id         = f2["neg_id"][()] #initial states of particles that develop negative p.a
nan_id         = f2["nan_id"][()] #initial states of particles that develop nan
high_id        = f2["high_id"][()] #initial states of particles that develop high p.a


print(neg_id)
print(high_id)
print(nan_id)


############################################### SAVE CSV FILES ######################################################33
#Initials to CSV file
header = ['id', 'aeq0_deg', 'lamda0_deg', 'ppar0', 'pper0', 'alpha0_deg', 'aeq00_deg','lamda00_deg','ppar00','pper00','alpha00_deg']
data = []
id=0
for aq0,l0,par0,per0,a0,aq00,l00,par00,per00,a00 in zip(aeq0, lamda0, ppar0, pper0, alpha0, aeq00, lamda00, ppar00, pper00, alpha00):
    data.extend([[id,aq0*R2D,l0*R2D,par0,per0,a0*R2D, aq00*R2D,l00*R2D,par00,per00,a00*R2D]])
    id=id+1 #id is the element's location in the list since the particles were saved like that
with open("dstr_data.csv", "w") as file1:
    writer = csv.writer(file1)
    writer.writerow(header)
    writer.writerows(data)
######################################## PLOT INITIAL DISTRIBUTION #######################################
#"""
fig, ax = plt.subplots()
for sec in range(0,sectors):
    ax.scatter(sec,aeq0_bins[sec],s=2,alpha=1)
ax.grid(alpha=.3)
ax.set(xlabel="Sectors",xlim=(0,sectors-1),xticks=np.arange(0,sectors),ylabel="dN",title="Aeq0 distribution, sector range "+str(sector_range)+" degrees")
ax.set_yscale("log")
fig.savefig("simulation_MM/"+str(population)+"_normals_1.png",dpi=200)

fig, ax = plt.subplots()
ax.scatter(lamda0*R2D,aeq0*R2D,s=0.5,alpha=0.1)
ax.grid(alpha=.3)
ax.set(xlabel="Latitude(deg)",ylabel="Equatorial P.A",title="Initial lat-aeq of simulated particles",ylim=(1,179),xlim=(-90,90),xticks=np.linspace(-90,90,5))
ax.axhline(y = 90, color ="b", linestyle="dashed")
fig.savefig("simulation_MM/"+str(sector_range)+"_normals_2.png",dpi=200)
#"""
