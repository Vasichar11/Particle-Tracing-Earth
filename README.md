# Project Description
This repository contains the code developed for my thesis titled 'Simulation of Wave-Particle Interaction using Parallel Processing'. This work was carried out during my studies in the Electrical and Computer Engineering Department at the Democritus University of Thrace.

To address the challenges of interacting with waves and handling a large population of particles, the following features were implemented:

- Different methods for generating particle distributions
- Parallel execution of the wave-particle interaction simulation in C++ using powerfull systems and exporting the results in hierarchical data format
- Loading these results in Python for further processing and analysis
- Extended visualization of the particle distributions using matplotlib library

In the end, the effect that waves with different intensities and frequencies have on different particle distributions during their propagation could be evindent.


## Table of Contents

- [Installation](#installation)
- [Workflow](#workflow)
- [Usage](#usage)
- [Acknowledgements](#acknowledgements)

## Workflow

To better understand the code structure, a flow diagram has been created. You can download the following PDF documents to access their hyperlinks and explore the corresponding GitHub modules:

- code.pdf
- Particle_distribution.pdf"

## Installation
There's support only for Linux at the moment. The intallation procedure is given below, step by step. Start with the dependencies, continue with cloning the repository, then install the python requirements and finally continue to [Usage](#usage).
### Dependencies
The commands are given for Debian, Ubuntu, and related distributions. If you have different distribution and/or different package manager, the commands are different but should be as simple as the ones that are listed below.
- libhdf5 and HighFive library (which is a git submodule to the repository, so you don't need to install explicitly) : ```sudo apt-get install libhdf5-dev```
- g++ compiler - with OpenMP support (if you want to run in parallel). g++ is most probably installed in your linux system  and ```echo |cpp -fopenmp -dM |grep -i open``` shows the version of OpenMP. It should give an output like "_OPENMP yyyymm", where yyyymm are the designations for the year and month of the OpenMP API version
- python3 - with PyQt5 support (if you want a GUI): python is most probably already installed, ```sudo apt-get update``` and ```sudo apt-get install python3``` and for PyQt5 ```pip install pyqt5```
### Clone Repository
Close the repository: ```git clone https://github.com/Vasichar11/Particle-Tracing-Earth```
Take not that the repository includes a git submodule for hdf5 file reading and writting operations:
https://github.com/BlueBrain/HighFive
### Requirements
There are a few requirements for python, listed in the requirements.txt file
- Create virtualenv: ```python3 -m venv my_venv``` 
- Activate: ```source my_venv/bin/activate``` 
- Install requirements (Particle-Tracing-Earth/requirements.txt): ```python3 -m pip install requirements.txt`` 


## Usage

### GUI
```python3 Particle-Tracing-Earth\tracer.py```
### CLI
1) ```make allclean``` cleans all simulation output data and object file
2) **Modify the constants.h header file** that includes all the parameter values for the simulation
3) ```make dstr <arg1> <arg2> <arg3> <arg4>``` will create the desired particle distribution
- ```<arg1>``` refers to the Pitch Angle (P.A) distribution 
- ```<arg2>``` refers to the Latitude distribution 
- ```<arg3>``` refers to the eta distribution 
- ```<arg4>``` refers to the Energy distribution 

These arguments take values from this list: (normal,uniform,evenly,constant). 
As the names suggest:

- normal is to create a normal distribution for the particles on the corresponding variable
- uniform is to create a uniform distribution for the particles on the corresponding variable
- evenly is to distribute the particles in evenly i.e. with a fixed step within a range on the corresponding variable
- constant is to distribute the particles in a constant value
4) ```make ray <stepsize>``` interpolates the ray values that are coming after ray tracing of the wave using the stepsize of the simulation. The stepsize is the timestep of the particles in the simulation. In every step we are checking if the particle is within the wave range so we need the values of the wave in these time frames.
5) ```make tracer <noWPI_time> <WPI_time>``` which is the WPI simulation command:
- <noWPI_time> is the time that the particles oscillate in adiavatic motion (without interacting with a wave). This time should be enough for the particles to be in a randomized state
- <WPI_time> is the time that the particles have to interact with the wave while bouncing in the Earth's magnetic field. Take note that this time should be much smaller than the noWPI_time sine the execution of the WPI code consumes more resources than the adiavatic motion simulation.
6) Post processing and visualization with python
- ```python3 src/visualization/Post_processing_and_Plots.py``` to make the binning process and visualize the results
- ```python3 src/visualization/distribution_plot.py``` to visualize the initial distribution
- ```python3 src/visualization/Ray_plot.py``` to visualize the ray electric and magnetic fields and latitude 
7) Visualization files have been created in the output/plots directory

## Acknowledgements
The initial simulation for single particle WPI was implemented serially in Python by PhD candidate Stelios Tourgaidis.



