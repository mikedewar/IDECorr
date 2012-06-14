#Author Parham Aram
#Date 13-06-2012
"""
This module provides the full neural field class, which describes a non-linear
integro-difference equation model and methods to simulate
the model
"""

#Standard library imports
from __future__ import division
import pylab as pb
import numpy as np
import time
import os
from scipy import signal
import scipy as sp

#My modules
from Bases2D import basis


class NF():

	"""class defining full neural model.

	Arguments
	----------
	kernel: 
			Connectivity kernel.

	sensor_kernel:

			output kernel, governs the sensor pick-up geometry
		
	obs_locns: List of matrix
				Sensors locations

	gamma : Gaussian basis function
				covariance function of the field disturbance
	gamma_weight: float
				amplitude of the field disturbance covariance function

	Sigma_varepsilon: matrix
				Observation noise covariance



	field_space: ndarray
				Spatiel field

	spacestep: float
			Spatial step size for descritization 	
	Attributes
	----------
	ny: int
		Number of sensors


	gen_ssmodel:
		Generate the full neural field model

	simulate:
		Simulate the full neural field model
	"""

        def __init__(self,kernel,sensor_kernel,obs_locns,gamma,gamma_weight,Sigma_varepsilon,simulation_space_x_y,spacestep):

                self.kernel = kernel
                self.sensor_kernel=sensor_kernel
                self.obs_locns=obs_locns
                self.gamma=gamma
                self.gamma_weight=gamma_weight
                self.Sigma_varepsilon=Sigma_varepsilon
                self.spacestep=spacestep
                self.simulation_space_x_y=simulation_space_x_y
                self.ny=len(self.obs_locns)
		

        def gen_ssmodel(self):

		"""
		generates full neural model

		Attributes:
		----------
		K: matrix
			matrix of connectivity kernel evaluated over the spatial domain of the kernel

		Sigma_e: matrix
			field disturbance covariance matrix
		Sigma_e_c: matrix
			cholesky decomposiotion of field disturbance covariance matrix
		Sigma_varepsilon_c: matrix
			cholesky decomposiotion of observation noise covariance matrix
		C: matrix
			matrix of sensors evaluated at each spatial location, it's not the same as C in the IDE model	

		"""
		print "generating full neural model"

		#Generate field meshgrid
                simulation_field_space_x,simulation_field_space_y=pb.meshgrid(self.simulation_space_x_y,self.simulation_space_x_y)


                K=0
                for i in range(len(self.kernel.Psi)):
                        K+=self.kernel.weights[i]*self.kernel.Psi[i](simulation_field_space_x,simulation_field_space_y)
        
                self.K=K


		#calculate field disturbance covariance matrix and its Cholesky decomposition

                gamma_space=pb.array(zip(simulation_field_space_x.flatten(),simulation_field_space_y.flatten()))

                N1,D1 = gamma_space.shape
                diff = gamma_space.reshape(N1,1,D1) - gamma_space.reshape(1,N1,D1)
                Sigma_e_temp=self.gamma_weight*np.exp(-np.sum(np.square(diff),-1)*(1./self.gamma.width))
                self.Sigma_e=Sigma_e_temp

                if hasattr(self,'Sigma_e_c'):
                        pass
                else:
                        self.Sigma_e_c=sp.linalg.cholesky(self.Sigma_e,lower=1)    

        #calculate Cholesky decomposition of observation noise covariance matrix
                Sigma_varepsilon_c=sp.linalg.cholesky(self.Sigma_varepsilon,lower=1)
                self.Sigma_varepsilon_c=Sigma_varepsilon_c

        #Calculate sensors at each spatial locations, it's not the same as C in the IDE model	
                t0=time.time()

                sensor_space=self.obs_locns
                N2,D2 = sensor_space.shape
                diff = sensor_space.reshape(N2,1,D2) - gamma_space.reshape(1,N1,D1)
                C=np.exp(-np.sum(np.square(diff),-1)*(1./self.sensor_kernel.width))
                self.C=C





	def simulate(self,T):

		"""Simulates the full neural field model

		Arguments
		----------

		T: ndarray
				simulation time instants
		Returns
		----------
		V: list of matrix
			each matrix is the neural field at a time instant

		Y: list of matrix
			each matrix is the observation vector corrupted with noise at a time instant
		"""

		Y=[]
		V=[]  

		spatial_location_num=(len(self.simulation_space_x_y))**2
		sim_field_space_len=len(self.simulation_space_x_y) 

		#initial field
		v0=pb.dot(self.Sigma_e_c,np.random.randn(spatial_location_num,1))
		v_membrane=pb.reshape(v0,(sim_field_space_len,sim_field_space_len))
        
		for t in T[1:]:
            
			v = pb.dot(self.Sigma_varepsilon_c,np.random.randn(len(self.obs_locns),1))
			w = pb.reshape(pb.dot(self.Sigma_e_c,np.random.randn(spatial_location_num,1)),(sim_field_space_len,sim_field_space_len))
			#print "simulation at time",t
			g=signal.fftconvolve(self.K,v_membrane,mode='same') 
			g*=(self.spacestep**2)
			v_membrane=g+w
            
			#Observation
			Y.append((self.spacestep**2)*(pb.dot(self.C,pb.reshape(v_membrane,(sim_field_space_len**2,1))))+v)
			V.append(v_membrane)

		return V,Y

	

