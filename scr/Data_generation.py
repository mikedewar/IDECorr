#Author Parham Aram
#Date 28-06-2010
'''This module simulates the full neural model and uses the data
in the state space model to estimate the states, 
the connectivity kernel parameters and the synaptic dynamics'''

#Standard library imports
from __future__ import division
from scipy import io
import pylab as pb
import numpy as np
import scipy as sp
from scipy import io
#My modules
from Bases2D import *
from NF import *
import IDEComponents
from circle import circle
from ObservationLocation import observation_locns
#space properties

dimension=2
#spatial step size
simulation_spacestep=0.5 # the step size should in a way that we have (0,0) in our kernel as the center
#field width
simulation_field_width=20
#spatial range
simulation_space_x_y=pb.arange(-(simulation_field_width)/2.,(simulation_field_width)/2.+simulation_spacestep,simulation_spacestep)
simulation_field_space_x,simulation_field_space_y=pb.meshgrid(simulation_space_x_y,simulation_space_x_y)
estimation_field_width=simulation_field_width
#time properties

#sampling rate  
Fs = 1e3                                       
#sampling period
Ts = 1/Fs   
t_end= 2.5
NSamples = t_end*Fs;
T = pb.linspace(0,t_end,NSamples);

#observation and field basis function locations

#-----------------------------------------------------------------------------------------------------------------------------------------------
#Define observation locations and basis function locations
Delta_s = 1.5#mm distance between sensors in mm
observation_locs_mm =observation_locns(2*simulation_spacestep,simulation_field_width,Delta_s)
#observation_locs_mm=pb.linspace(-8.,8.,14)
print 'Sensor centers:',observation_locs_mm
obs_locns_x,obs_locns_y=pb.meshgrid(observation_locs_mm,observation_locs_mm)
obs_locns=pb.array(zip(obs_locns_x.flatten(),obs_locns_y.flatten()))
#Define connectivity kernel basis functions

#centers
psi0_center=pb.matrix([[-0.5],[-0.5]]);psi1_center=pb.matrix([[0],[0]]);psi2_center=pb.matrix([[0.5],[0.5]])
#widths
#psi0_width=1.8**2;psi1_width=2.4**2;psi2_width=6.**2
psi0_width=2.4**2;psi1_width=2.4**2;psi2_width=2.4**2

#weights
psi0_weight=-.07;psi1_weight=.07;psi2_weight=.02
#define each kernel basis functions
psi0=basis(psi0_center,psi0_width,dimension)
psi1=basis(psi1_center,psi1_width,dimension)
psi2=basis(psi2_center,psi2_width,dimension)

NF_Connectivity_kernel_basis_functions=[psi0,psi1,psi2]
#NF_Connectivity_kernel_basis_functions=[psi0]

NF_Connectivity_kernel_weights=pb.array([psi0_weight,psi1_weight,psi2_weight])
#NF_Connectivity_kernel_weights=pb.array([psi0_weight])

NF_Connectivity_kernel=IDEComponents.Kernel(NF_Connectivity_kernel_basis_functions,NF_Connectivity_kernel_weights)


#Field noise
gamma_center=pb.matrix([[0],[0]])
gamma_width=1.3**2 
gamma_weight=0.1
gamma=basis(gamma_center,gamma_width,dimension)
#Observation noise
varepsilon=0.1
Sigma_varepsilon =varepsilon*np.eye(len(obs_locns),len(obs_locns))

#Define sensor geometry
sensor_center=pb.matrix([[0],[0]])
sensor_width=0.9**2 
sensor_kernel=basis(sensor_center,sensor_width,dimension)
#plot sensors and basis functions
#l1=[circle(cent,2*pb.sqrt(pb.log(2),)*pb.sqrt(sensor_width)/2,'b') for cent in obs_locns]
#pb.show()

#ignore first 100 observations allowing the model's initial transients to die out
First_n_observations=100
number_of_iterations=50
#populate the model
NF_model=NF(NF_Connectivity_kernel,sensor_kernel,obs_locns,gamma,gamma_weight,Sigma_varepsilon,simulation_space_x_y,simulation_spacestep)
#generate the Neural Field model
NF_model.gen_ssmodel()
V,Y=NF_model.simulate(T)
data=dict(Y=Y)
io.savemat('data',data)

