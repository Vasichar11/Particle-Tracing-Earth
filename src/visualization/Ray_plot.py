import matplotlib.pyplot as plt
import numpy as np
import h5py
import os
import re

# Read header file to get pulse duration and pwr
with open('headers/constants.h', 'r') as file:
    header_code = file.read()
# Extract values
match1 = re.search(r"const\s+real\s+pulse_duration\s*=\s*(.*?);", header_code)
match2 = re.search(r"const\s+real\s+pwr\s*=\s*(.*?);", header_code)
if match1:
    pulse_duration = float(match1.group(1))
    print("Pulse duration:", pulse_duration)
else:
    raise ValueError("pulse_duration not found in constants.h")   

if match2:
    pwr_expr = match2.group(1)
    try:
        pwr = float(eval(pwr_expr))
        print("Power:", pwr)
    except (SyntaxError, ValueError) as e:
        print("Error evaluating power expression:", str(e))
else:
    raise ValueError("pwr not found in constants.h")


# Read particle distribution from file entered by user
files_dir = 'output/files'  # Path to the output directory
plots_dir = 'output/plots'  # Path to the output directory
h5_files = [file for file in os.listdir(files_dir) if (file.startswith('interpolated_ray') and file.endswith('.h5'))]

# Display the list of .h5 files
print("Available ray .h5 files:")
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
f1 = h5py.File(selected_file,"r")
lat = f1["lat_int"][()]
Ezw   = f1["Ezw"][()]  
Bw_ray = f1["Bw_ray"][()]
Bzw = f1["Bzw"][()]
Ezw = f1["Ezw"][()]
time = f1["time"][()]

fig, ax = plt.subplots(2, 2)

ax[0,0].plot(time, Bw_ray, label="Bw_ray")
ax[1,1].set(xlabel="time(s)",ylabel="Field(T)",ylim=(min(Bw_ray),max(Bw_ray)))
ax[0,0].legend()

ax[0,1].plot(time, Ezw, label="Ezw")
ax[1,1].set(xlabel="time(s)",ylabel="Field(V/m)",ylim=(min(Ezw),max(Ezw)))
ax[0,1].legend()

ax[1,0].plot(time, Bzw, label="Bzw")
ax[1,1].set(xlabel="time(s)",ylabel="Field(T)",ylim=(min(Bzw),max(Bzw)))
ax[1,0].legend()

ax[1,1].plot(time, lat, label="lat")
ax[1,1].set(xlabel="time(s)",ylabel="Field(deg)",ylim=(min(lat),max(lat)))
ax[1,1].legend()

fig.savefig(os.path.join(plots_dir,selected_file.split("/")[-1].split(".h5")[0]+".png"))