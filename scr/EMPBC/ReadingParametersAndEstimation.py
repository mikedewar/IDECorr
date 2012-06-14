#Built-in modules
from __future__ import division
import pylab as pb
import numpy as np
import scipy as sp
from scipy import io
import os
#My modules
from Bases2D import *
from IDEModel import *
import IDEComponents
from circle import circle
import ActivationFunction
#simulation properties
#----------------------------------------------------------------------------------------------------------------------------------------------
dimension=2
#step size
Spacestep=0.1 # mm
#field width
FieldWidth=3.6 # the width of the field in mm
#spatial range
estimation_space_x_y=pb.arange(-(FieldWidth)/2.,(FieldWidth)/2.+Spacestep,Spacestep)
# Define Sampling properties
Fs = 5e3   #sampling rate   Hz                                    
Ts = 1/Fs   #sampling period, second
#t_end= 0.4 # seconds
#NSamples = t_end*Fs;
#T = pb.linspace(0,t_end,NSamples);
#-----------------------------------------------------------------------------------------------------------------------------------------------
#Define observation locations and basis function locations
NumberOfSensors_x_y=10
observation_locs_mm=pb.linspace(-FieldWidth/2,FieldWidth/2,NumberOfSensors_x_y)
print 'Sensor centers:',observation_locs_mm
obs_locns_x,obs_locns_y=pb.meshgrid(observation_locs_mm,observation_locs_mm)
obs_locns=pb.array(zip(obs_locns_x.ravel(),pb.flipud(obs_locns_y).ravel()))
#remove unwanted channels
NoisyChannelsIndex=[0,9,12, 14, 18, 23, 25, 34, 41, 45, 46, 53, 55, 56, 66, 80, 84,86,90,99]
ChannelIndex=range(100)
CleanChannelIndex=pb.setxor1d(NoisyChannelsIndex,ChannelIndex)
#Observation locations are arranged in a way that it starts from right top corner, and moves horizontally
obs_locns=obs_locns[CleanChannelIndex,:]
sensor_FWHM=0.6#Full width at half maximum mm #between 400 and 800 um
#[circle(i,sensor_FWHM/2.,'b') for i in obs_locns]
#Define sensor geometry
sensor_center=pb.matrix([[0],[0]])
sensor_width=float((sensor_FWHM/(2*pb.sqrt(pb.log(2))))**2) #mm^2
sensor_kernel=basis(sensor_center,sensor_width,dimension)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Define basis functions locations and widths
NBasisFunction_x_y=8
S=pb.linspace(-FieldWidth/2.,FieldWidth/2.,NBasisFunction_x_y)
phi_centers_x,phi_centers_y=pb.meshgrid(S,S)
#basis functions are arranged in a way that it starts from right top corner, and moves horizontally
phi_centers=zip(phi_centers_x.flatten(),pb.flipud(phi_centers_y).flatten())

#Define basis functions
phi_FWHM=0.85  #must be float mm
phi_width=float((phi_FWHM/(2*pb.sqrt(pb.log(2))))**2) #mm^2
#place field basis functions in an array in the form of n_x*1
Phi=pb.array([basis(center,phi_width,dimension) for center in phi_centers],ndmin=2).T
field=IDEComponents.Field(Phi)
#[circle(i,phi_FWHM/2.,'k') for i in phi_centers]
#------------------------------------------------------------------------------------------------------------------------------------------------
#Define connectivity kernel basis functions
#Connectivity kernel basis functions' centers
psi0_center=pb.matrix([[0],[0]])
psi1_center=pb.matrix([[0],[0]])
#Kernel basis functions' widths
psi0_FWHM=1.2 #Full width at half maximum mm
psi1_FWHM=1.8 #Full width at half maximum mm
psi0_width=float((psi0_FWHM/(2*pb.sqrt(pb.log(2))))**2) #mm^2
psi1_width=float((psi1_FWHM/(2*pb.sqrt(pb.log(2))))**2) #mm^2
#Kernel basis functions' widths
psi0_weight=0
psi1_weight=0
psi0=basis(psi0_center,psi0_width,dimension)
psi1=basis(psi1_center,psi1_width,dimension)
#create a list of connectivity kernel basis functions
Connectivity_kernel_basis_functions=[psi0,psi1]
Connectivity_kernel_weights=pb.array([psi0_weight,psi1_weight])
Connectivity_kernel=IDEComponents.Kernel(Connectivity_kernel_basis_functions,Connectivity_kernel_weights)
#Define field covariance function and observation noise
#Field noise
gamma_center=pb.matrix([[0],[0]])
gamma_width=.7**2 #1.1
gamma_weight=.2
gamma=basis(gamma_center,gamma_width,dimension)
#Observation noise
varepsilon=.2 #10
Sigma_varepsilon =varepsilon*np.eye(len(obs_locns),len(obs_locns))

#Define sigmoidal activation function and zeta
#------------------------------------
fmax=1
v0=1.8
varsigma=0.56
act_fun=ActivationFunction.ActivationFunction(v0,fmax,varsigma)
zeta=100
#----------Field initialasation and number of iterations for estimation----------------------------
mean=[0]*len(Phi)
P0=10*pb.eye(len(mean))
x0=pb.multivariate_normal(mean,P0,[1]).T
number_of_iterations=50
model=IDE(Connectivity_kernel,field,sensor_kernel,obs_locns,gamma,gamma_weight,Sigma_varepsilon,act_fun,x0,P0,zeta,Ts,estimation_space_x_y,Spacestep)
model.gen_ssmodel()
Y_temp=io.loadmat('Y0')
Y=Y_temp['Data']
ps_estimate=para_state_estimation(model,order=2)
ps_estimate.itrerative_state_parameter_estimation(Y,number_of_iterations)


