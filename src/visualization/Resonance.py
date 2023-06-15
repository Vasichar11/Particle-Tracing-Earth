import numpy as np 
import sys
import matplotlib.pyplot as plt
import h5py
import math
import matplotlib.animation as manimation
from random import randint
import threading
import os


D2R=np.pi/180
R2D=1/D2R
f1 = h5py.File("output/files/43200p_wpi.h5","r")
saved_id = f1["saved_id"][()]
saved_max_dEkin = f1["saved_max_dEkin"][()]
saved_maxEkin_time = f1["saved_maxEkin_time"][()]
saved_max_dPA = f1["saved_max_dPA"][()]
saved_maxdPA_time = f1["saved_maxdPA_time"][()]
f1.close()



