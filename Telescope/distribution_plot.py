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
f1 = h5py.File("h5files/1000p_normalAEQ_normalLAMDA.h5","r")
aeq0         = f1["aeq"][()]
lamda0       = f1["lat"][()]
aeq0_bins    = f1["aeq0_bins"][()]
f1.close()

f2 = h5py.File("h5files/1000p_no_wpi.h5","r")
lamda00        = f2["lamda00"][()] #states when noWPI stops
ppar00         = f2["ppar00"][()]
pper00         = f2["pper00"][()]
alpha00        = f2["alpha00"][()]
aeq00          = f2["aeq00"][()]
eta00          = f2["eta00"][()]
time00         = f2["time00"][()]


#Save initials to CSV file
header = ['id', 'aeq0_deg', 'lamda0_deg','aeq00_deg','lamda00_deg']
data = []
p=0
for a0,l0,a00,l00 in zip(aeq0, lamda0, alpha00, lamda00):
    data.extend([[p,a0*R2D,l0*R2D,a00*R2D,l00*R2D]])
    p=p+1 #id is the element's location in the list since the particles were saved like that
with open("dstr_before_and_after_WPI.csv", "w") as f:
    writer = csv.writer(f)
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
fig.savefig("simulation_MM/10e6_normals.png",dpi=200)

fig, ax = plt.subplots()
ax.scatter(lamda0*R2D,aeq0*R2D,s=0.5,alpha=0.1)
ax.grid(alpha=.3)
ax.set(xlabel="Latitude(deg)",ylabel="Equatorial P.A",title="Initial lat-aeq of simulated particles",ylim=(1,179),xlim=(-90,90),xticks=np.linspace(-90,90,5))
ax.axhline(y = 90, color ="b", linestyle="dashed")
fig.savefig("simulation_MM/10e6_normals_2.png",dpi=200)
#"""