#Built-in modules
from __future__ import division
import pylab as pb
import numpy as np
import time
import os
from scipy import signal
import scipy as sp
import UKF
import LS
class IDE():

	"""class defining a non-linear, Integro-Difference, discrete-time state space model.

	Arguments
	----------
	kernel: IDEComponents.Kernel instance
			IDE kernel.
	field : IDEComponents.Field instance
			spatial field.
	sensor_kernel: Bases2D.basis instance
			output kernel, governs the sensor pick-up geometry.	
	observation_locs_mm:ndarray
				Sensors locations along x or y axes

	gamma : Bases2D.basis instance
				covariance function of the field disturbance
	gamma_weight: float
				amplitude of the field disturbance covariance function

	Sigma_varepsilon: ndarry
				Observation noise covariance
	act_fun :ActivationFunction.ActivationFunction instance
				Sigmoidal activation function
	x0: ndarray
				Initial state
	P0: ndarray
				Initial state covariance matrix
	zeta: float
				Inverse synaptic time constant
	Ts: float
				Sampling time

	estimation_space_x_y: ndarray
				Spatiel field

	spacestep: float
			Spatial step size for descritization 
	
	Attributes
	----------
	nx: int
		Number of basis functions
	ny: int
		Number of sensors

	xi: float
		Time constant parameter

	gen_ssmodel:
		Generate the state space model

	simulate:
		Simulate the state space model

	state_equation:
		Evaluate state transition function at a given state

	"""


	def __init__(self,kernel,field,sensor_kernel,obs_locns,gamma,gamma_weight,Sigma_varepsilon,act_fun,x0,P0,zeta,Ts,estimation_space_x_y,spacestep):

		self.kernel = kernel
		self.field = field
		self.sensor_kernel=sensor_kernel
		self.obs_locns=obs_locns
		self.gamma=gamma
		self.gamma_weight=gamma_weight
		self.Sigma_varepsilon=Sigma_varepsilon
		self.act_fun = act_fun
		self.x0=x0
		self.P0=P0
		self.zeta=zeta
		self.Ts=Ts
		self.spacestep=spacestep
		self.estimation_space_x_y=estimation_space_x_y
		self.nx=len(self.field.Phi)
		self.ny=len(obs_locns)
		self.n_theta=len(self.kernel.Psi)
		self.xi=1-(self.Ts*self.zeta)

	def gen_ssmodel(self):

		'''Generates non-linear, Integro-Difference, discrete-time state space model.

		Atributes:
		----------
		Gamma: ndarray
			Inner product of field basis functions

		Gamma_inv: ndarray
			Inverse of Gamma

		Sigma_e: ndarray
			Covariance matrix of the field disturbance

		Sigma_e_inv: ndarray
			Inverse of Sigma_e

		Sigma_e_c: ndarray
			cholesky decomposiotion of field disturbance covariance matrix
		
		Sigma_varepsilon_c: ndarray
			cholesky decomposiotion of observation noise covariance matrix

		Phi_values: ndarray  nx by number of spatial locations
			field basis functions values over the spatial field

		psi_conv_Phi_values: dictioanary of ndarray ; each ndarray is nx by number of spatial locations
			convolution of the  connectivity kernel basis functions with the field basis functions
			evaluate over space


		Gamma_inv_psi_conv_Phi:dictionary of ndarray; each ndarray is nx by number of spatial locations
			the product of inverse gamme withe psi0_conv_Phi


		lambda_i:dictionary of list of ndarray; each list with the length of nx and each ndarray is number of spatial locations by nx


		lambda_tilde_i:dictionary of list of ndarray; each list with the length of nx and each ndarray is nx by number of spatial locations 


		C: ndarray
			Observation matrix

		Gamma_inv_Psi_conv_Phi: ndarray nx by number of spatial locations
			the convolution of the kernel  with field basis functions at discritized spatial points

		'''

		#Generate Gamma
		if hasattr(self,'Gamma'):
			pass
		else:
			t_total=time.time()
			t_Gamma=time.time()
			#calculate Gamma=PhixPhi.T; inner product of the field basis functions
			Gamma=pb.dot(self.field.Phi,self.field.Phi.T)
			Gamma_inv=pb.inv(Gamma)
			self.Gamma=Gamma.astype('float')
			self.Gamma_inv=Gamma_inv
			print 'Elapsed time in seconds to calculate Gamma inverse is',time.time()-t_Gamma

		#Generate field covariance matrix
		if hasattr(self,'Sigma_e'):
			pass
		else:
			print "calculating field noise covariances"
			t_Sigma_e_c=time.time()
			gamma_convolution_vecrorized=pb.vectorize(self.gamma.conv)
			gamma_conv_Phi=gamma_convolution_vecrorized(self.field.Phi).T 
			#[gamma*phi1 gamma*phi2 ... gamma*phin] 1 by nx
			Pi=pb.dot(self.field.Phi,gamma_conv_Phi) #nx by nx ndarray
			Pi=Pi.astype('float')
			Sigma_e=pb.dot(pb.dot(self.gamma_weight*Gamma_inv,Pi),Gamma_inv.T)
			Sigma_e_c=sp.linalg.cholesky(Sigma_e,lower=1)
			self.Sigma_e=Sigma_e
			self.Sigma_e_inv=pb.inv(Sigma_e)
			self.Sigma_e_c=Sigma_e_c
			Sigma_varepsilon_c=sp.linalg.cholesky(self.Sigma_varepsilon,lower=1)
			self.Sigma_varepsilon_c=Sigma_varepsilon_c

			print 'Elapsed time in seconds to calculate Cholesky decomposition of the noise covariance matrix',time.time()-t_Sigma_e_c


		if hasattr(self,'Phi_values'):
			pass
		else:

		
			#Generate field meshgrid
			estimation_field_space_x,estimation_field_space_y=pb.meshgrid(self.estimation_space_x_y,self.estimation_space_x_y)
			estimation_field_space_x=estimation_field_space_x.ravel()
			estimation_field_space_y=estimation_field_space_y.ravel()
			#calculate Phi_values
			#Phi_values is like [[phi1(r1) phi1(r2)...phi1(rn_r)],[phi2(r1) phi2(r2)...phi2(rn_r)],..[phiL(r1) phiL(r2)...phiL(rn_r)]]
			Phi_values=[self.field.Phi[i,0](estimation_field_space_x,estimation_field_space_y) for i in range(self.nx)]
			self.Phi_values=pb.squeeze(Phi_values) #it's nx by number of apatial points

			#vectorizing kernel convolution method
			psi_convolution_vectorized=pb.empty((self.n_theta,1),dtype=object)
			for i in range(self.n_theta):
				psi_convolution_vectorized[i,0]=pb.vectorize(self.kernel.Psi[i].conv)

			#find convolution between kernel and field basis functions analytically
			psi_conv_Phi=pb.empty((self.nx,self.n_theta),dtype=object)#nx by n_theta
			for i in range(self.n_theta):
				psi_conv_Phi[:,i]=psi_convolution_vectorized[i,0](self.field.Phi).ravel()

			self.psi_conv_Phi_analytic=psi_conv_Phi

			#ecaluate convolution between kernel and field basis functions at spatial locations
			psi_conv_Phi_values=pb.empty((self.n_theta,self.nx,len(self.estimation_space_x_y)**2),dtype=float)
			for i in range(self.n_theta):
				for j in range(self.nx):
					psi_conv_Phi_values[i,j,:]=psi_conv_Phi[j,i](estimation_field_space_x,estimation_field_space_y)
			self.psi_conv_Phi_values=psi_conv_Phi_values


			#It is nx by n_theta by n_r psi0->[:,0,:] psi1->[:,0,:] ... psi_ntheta->[:,ntheta,:]
			#(checked)
			self.Gamma_inv_psi_conv_Phi=pb.dot(self.Gamma_inv,self.psi_conv_Phi_values)

			#find Sigma_{e-1}*Gamma^{-1}
			Sigma_e_inv_Gamma_inv_T=pb.dot(self.Sigma_e_inv,self.Gamma_inv)
			#find bits in second term of Xi0 and Xi3 : nx x n_theta x nr
			#(checked)
			SIGI_temp=pb.dot(Sigma_e_inv_Gamma_inv_T,self.psi_conv_Phi_values)
			#I transpose it because I need it's transpose in parameter estimation
			self.SIGI=pb.transpose(SIGI_temp,(1,0,2))
			#Calculate bits in lambda_i lambda_tilde_i
			#find Gamma^{-T}*Sigma_{e}^{-1}
			Gamma_inv_T_Sigma_e_inv=pb.dot(self.Gamma_inv.T,self.Sigma_e_inv)



			product_phi_i_psi=pb.empty((self.nx,self.nx,self.n_theta),dtype=object)#nx by nx by ntheta
			product_phi_i_psi_values=pb.empty((self.nx,self.nx,self.n_theta,len(self.estimation_space_x_y)**2),dtype=float)	
			#vectorizing field basis functions' product method for speed, we override __or__ to perform Gaussian basis functions product
			for i in range(self.nx):
				product_vectorized=pb.vectorize(self.field.Phi[i,0].__or__)	
				product_phi_i_psi[i,:,:]=product_vectorized(psi_conv_Phi)
				for j in range(self.nx):
					for k in range(self.n_theta):
						product_phi_i_psi_values[i,j,k,:]=product_phi_i_psi[i,j,k](estimation_field_space_x,estimation_field_space_y)


			product_phi_i_psi_values_temp=pb.transpose(product_phi_i_psi_values,(0,2,3,1))
			t0=time.time()
			lambda_i=pb.dot(product_phi_i_psi_values_temp,Gamma_inv_T_Sigma_e_inv)
			#I transpose it for the speed in parameter estimation, if I use transpose it's twice faster in parameter estimation
			#it should be nx x n_theta x nr x nx but it is transpoded as nx x nx x nr x n_theta
			# the result in parameter estimation part must be transposed
			self.lambda_i=pb.transpose(lambda_i,(0,3,2,1))
			print 'lambda_i',time.time()-t0
			
			t0=time.time()
			lambda_tilde_i=pb.dot(product_phi_i_psi_values_temp,self.Gamma_inv.T)
			self.lambda_tilde_i=pb.transpose(lambda_tilde_i,(0,3,2,1))
			print 'lambda_tilda_i',time.time()-t0
			#calculate bits in the second term of Xi0 and Xi3
					

		


			#Calculate observation matrix
			
		if hasattr(self,'C'):
			pass
		else:
			#Generate Observation locations grid
			obs_locns_x=self.obs_locns[:,0]
			obs_locns_y=self.obs_locns[:,1]
			t_observation_matrix=time.time()
			sensor_kernel_convolution_vecrorized=pb.vectorize(self.sensor_kernel.conv)
			sensor_kernel_conv_Phi=sensor_kernel_convolution_vecrorized(self.field.Phi).T #first row 
			#[m*phi_1 m*phi2 ... m*phin]
			C=pb.empty(([self.ny,self.nx]),dtype=float)
			for i in range(self.nx):
				C[:,i]=sensor_kernel_conv_Phi[0,i](obs_locns_x,obs_locns_y)

			self.C=C
			print 'Elapsed time in seconds to calculate observation matrix C is',time.time()-t_observation_matrix	
			print 'Elapsed time in seconds to generate the model',time.time()-t_total


		#We need to calculate this bit at each iteration	
		##Finding the convolution of the kernel  with field basis functions at discritized spatial points
		self.Gamma_inv_Psi_conv_Phi=pb.sum(self.kernel.weights[pb.newaxis,:,pb.newaxis]*self.Gamma_inv_psi_conv_Phi,axis=1)



	def simulate(self,T):

		"""
		generates nonlinear IDE


		Arguments
		----------
		init_field: ndarray
				initial state in a form of nx x 1
		T: ndarray
				simulation time instants
		Returns
		----------
		X: list of ndarray
			each ndarray is the state vector at a time instant

		Y: list of ndarray
			each ndarray is the observation vector corrupted with noise at a time instant
		"""


		Y = []		#at t=0 we don't have observation
		X = []		#I don't save the initial state and X starts at t=1
		x=self.x0
		firing_rate_temp=pb.dot(x.T,self.Phi_values)
 		firing_rate=self.act_fun.fmax/(1.+pb.exp(self.act_fun.varsigma*(self.act_fun.v0-firing_rate_temp))).T		
		print "iterating"
		for t in T[1:]:
			w = pb.dot(self.Sigma_e_c,np.random.randn(self.nx,1))
			v = pb.dot(self.Sigma_varepsilon_c,np.random.randn(len(self.obs_locns),1))
			print "simulation at time",t
			g=pb.dot(self.Gamma_inv_Psi_conv_Phi,firing_rate)
			g *=(self.spacestep**2)
			x=self.Ts*g+self.xi*x+w

			X.append(x)
			Y.append(pb.dot(self.C,x)+v)
			firing_rate_temp=pb.dot(x.T,self.Phi_values)
			firing_rate=self.act_fun.fmax/(1.+pb.exp(self.act_fun.varsigma*(self.act_fun.v0-firing_rate_temp))).T		

		return X,Y


	def  state_equation(self,x):

		'''state equation for sigma points propogation '''
		firing_rate_temp=pb.dot(x.T,self.Phi_values)
 		firing_rate=self.act_fun.fmax/(1.+pb.exp(self.act_fun.varsigma*(self.act_fun.v0-firing_rate_temp))).T	
		g=pb.dot(self.Gamma_inv_Psi_conv_Phi,firing_rate)
		g *=(self.spacestep**2)
		x=self.Ts*g+self.xi*x
		return x



class para_state_estimation():

	def __init__(self,model,order):

		'''this is to estimate state and connectivity kernel parameters

		Arguments:
		----------
			model: IDE instance
			order: 1 or 2
				specify the zero or first Taylor approximation to the non-linearity '''

		self.model=model
		self.order=order


	def estimate_kernel(self,X,Y,P,M):

		"""
			estimate the ide model's kernel weights using Least Square method
	
			Arguments
			----------
			X: list of matrix
				state vectors

			Returns
			---------
			Least Square estimation of the IDE parameters, see the corresponding pdf file
		"""



		if self.order==2:


			#Q, Lambda and Lambda_tilde are already multiplied by Ts

			#Calculate R1,R2 and R3

			T=len(X)
			

			R1=0;R2=0;R3=0;Xi0=0;Xi1=0;Xi2=0;Xi3=0;Xi4=0;Sigma_varepsilon=0
			#nx byn_theta by nr
			Psi=self.model.Gamma_inv_psi_conv_Phi
			#I transpose it here for speed, to multiply firing rate with Psi pointwise and add over space
			#nx by nr by n_theta
			Psi_T=pb.transpose(self.model.Gamma_inv_psi_conv_Phi,(0,2,1))

			for t in range(T-1):

				
				# 1 x n_r [v(r1),v(r2),...v(rn)]
				firing_rate_temp=pb.dot(X[t].T,self.model.Phi_values)
				# 1 x n_r [f(v(r1)),f(v(r2)),...f(v(rn))]
				firing_rate=self.model.act_fun.fmax/(1.+pb.exp(self.model.act_fun.varsigma*(self.model.act_fun.v0-firing_rate_temp)))	
				#1 x n_r [f'(v(r1)),f'(v(r2)),...f'(v(rn))]
				firing_rate_derivative=self.model.act_fun.varsigma*firing_rate*(1-(firing_rate/self.model.act_fun.fmax))

				#calculate q
				#1 by n_x by n_theta
				g=pb.dot(firing_rate,Psi_T)

				g *=(self.model.spacestep**2)	
				q=self.model.Ts*g
				#This is needed for the first term in Xi0, Xi2 and Xi3
				# it's 1 by nx by n_theta and I get rid of the first axis
				q=q.reshape(self.model.nx,self.model.n_theta)
				temp=pb.dot(self.model.Sigma_e_inv,q)#nx by ntheta



				#bits in the second term of Xi0 and Xi3
				#SIGI:ntheta by nx by nr
				#Xi_0_3:ntheta by nx by nr
				Xi_0_3=self.model.SIGI*firing_rate_derivative


				#calculate Lambda and Lambda_tilde
				#pointwise multiplication between each lambda components and derivative of firing rate for psi0,psi1 and psi2 and add them up

				Lambda=pb.dot(firing_rate_derivative,self.model.lambda_i)#1 by nx by nx by ntheta
				#here we transpose because the bit in model generation was transposed for speed
				Lambda=pb.transpose(Lambda,(0,1,3,2))*(self.model.spacestep**2)*self.model.Ts#1 by nx by ntheta by nx
				
				Lambda_tilde=pb.dot(firing_rate_derivative,self.model.lambda_tilde_i)*(self.model.spacestep**2)*self.model.Ts
				#1 by nx by nx by ntheta
				#calculate Xi_0
				#loop is faster than broadcasting
				Xi0_temp=pb.empty((1,self.model.n_theta),dtype=float)
				for i in range(self.model.n_theta):
					Xi0_temp[:,i]=pb.sum(pb.sum(pb.dot(M[t+1],Xi_0_3[i,:,:])*self.model.Phi_values,axis=1),axis=0)

				Xi0_temp=Xi0_temp*(self.model.spacestep**2)*self.model.Ts
				Xi0+=pb.dot(X[t+1].T,temp)+Xi0_temp

				#calculate Xi1

				Xi1+=M[t+1]+pb.dot(X[t],X[t+1].T)

				#calculate Xi2
				Lambda_dot_P_t=pb.dot(Lambda,P[t])
				Xi2_temp=[pb.dot(Lambda_dot_P_t[0,i],Lambda_tilde[0,i]) for i in range(len(Lambda_tilde))]
				Xi2+=pb.dot(q.T,temp)+pb.sum(Xi2_temp,axis=0)

				#calculate Xi3
				#loop is faster than broadcasting
				Xi3_temp=pb.empty((1,self.model.n_theta),dtype=float)
				for i in range(self.model.n_theta):
					Xi3_temp[:,i]=pb.sum(pb.sum(pb.dot(P[t],Xi_0_3[i,:,:])*self.model.Phi_values,axis=1),axis=0)

				Xi3_temp=Xi3_temp*(self.model.spacestep**2)*self.model.Ts
			
				Xi3+=pb.dot(X[t].T,temp)+Xi3_temp

				#calculate Xi4
				Xi4+=P[t]+pb.dot(X[t],X[t].T)
				#calculate observation noise covariance
				Sigma_varepsilon+=pb.dot(Y[t+1]-pb.dot(self.model.C,X[t+1]),(Y[t+1]-pb.dot(self.model.C,X[t+1])).T)+dots(self.model.C,P[t+1],self.model.C.T)

			#calculate r1,r2,r'1,r'2
			Sigma_varepsilon=Sigma_varepsilon/(T-1)
			r1=pb.trace(pb.dot(Xi1,self.model.Sigma_e_inv))
			r2=pb.trace(pb.dot(Xi4,self.model.Sigma_e_inv))
			rdash1=dots(Xi3,pb.inv(Xi2),Xi0.T)
			rdash2=dots(Xi3,pb.inv(Xi2),Xi3.T)
			xi=float((r1-rdash1)/(r2-rdash2))

			theta=pb.dot(pb.inv(Xi2),Xi0.T-xi*Xi3.T)
			parameters=[float(theta[i]) for i in range(theta.shape[0])]
			parameters.append(float(xi))
			return parameters,Sigma_varepsilon





	def itrerative_state_parameter_estimation(self,Y,max_it):

		"""Two part iterative algorithm, consisting of a state estimation step followed by a
		parameter estimation step
		
		Arguments:
		---------
		Y: list of matrices
			Observation vectors
		max_it: int
			maximum number of iterations """


		xi_est=[]
		kernel_weights_est=[]
		Sigmavar_epsilon_est=[]
		# generate a random state sequence
		Xb= [np.random.rand(self.model.nx,1) for t in Y]
		Pb=[pb.zeros((Xb[0].shape[0],Xb[0].shape[0]))]*len(Xb)
		Mb=[pb.zeros((Xb[0].shape[0],Xb[0].shape[0]))]*len(Xb)
		# iterate
		keep_going = 1
		it_count = 0
		print " Estimatiing IDE's kernel, the time constant parameter, observation noise covariance and field weights"
		t0=time.time()
		while keep_going:
			
			temp,Sigma_varepsilon=self.estimate_kernel(Xb,Y,Pb,Mb)
			xi_est.append(float(temp[-1]))
			kernel_weights_est.append(temp[0:-1])
			#Sigmavar_epsilon_est.append(Sigma_varepsilon)
			self.model.kernel.weights,self.model.xi=pb.array(temp[0:-1]),float(temp[-1])
			#self.model.Sigma_varepsilon=Sigma_varepsilon
			#self.model.Sigma_varepsilon_c=sp.linalg.cholesky(Sigma_varepsilon,lower=1)
			self.model.gen_ssmodel()
			ukf_instance=UKF.ukf(self.model.nx,self.model.x0,self.model.P0,self.model.C,self.model.Sigma_e,self.model.Sigma_varepsilon,self.model.state_equation,kappa=0.0,alpha_sigma_points=1e-3,beta_sigma_points=2)
			Xb,Pb,Mb,Xf,Pf=ukf_instance.rtssmooth(Y)
			self.model.Xb=Xb
			self.model.Pb=Pb
			self.model.Xf=Xf
			self.model.Pf=Pf
			self.model.xi_est=xi_est
			self.model.kernel_weights_est=kernel_weights_est
			self.model.Sigmavar_epsilon_est=Sigmavar_epsilon_est
			self.model.x0=Xb[0]
			self.model.P0=Pb[0]
			print it_count, " Kernel current estimate: ", self.model.kernel.weights, "xi", self.model.xi#,"Sigma_varepsilon",self.model.Sigma_varepsilon
			if it_count == max_it:
				keep_going = 0
			it_count += 1
		print "Elapsed time in seconds is", time.time()-t0


	def LS_itrerative_state_parameter_estimation(self,Y,max_it):

		"""Two part iterative algorithm, consisting of a state estimation step followed by a
		parameter estimation step
		
		Arguments:
		---------
		Y: list of matrices
			Observation vectors
		max_it: int
			maximum number of iterations """


		xi_est=[]
		kernel_weights_est=[]
		# generate a random state sequence
		Xb= [np.random.rand(self.model.nx,1) for t in Y]
		# iterate
		keep_going = 1
		it_count = 0
		print " Estimatiing IDE's kernel, the time constant parameter, observation noise covariance and field weights"
		t0=time.time()
		while keep_going:
			Ls_instance=LS. para_state_estimation(self.model)
			temp=Ls_instance.estimate_kernel(Xb)
			xi_est.append(float(temp[-1]))
			kernel_weights_est.append(temp[0:-1])
			self.model.kernel.weights,self.model.xi=pb.array(temp[0:-1]),float(temp[-1])
			self.model.gen_ssmodel()
			ukf_instance=UKF.ukf(self.model.nx,self.model.x0,self.model.P0,self.model.C,self.model.Sigma_e,self.model.Sigma_varepsilon,self.model.state_equation,kappa=0.0,alpha_sigma_points=1e-3,beta_sigma_points=2)
			Xb,Pb,Mb,Xf,Pf=ukf_instance.rtssmooth(Y)
			self.model.Xb=Xb
			self.model.Pb=Pb
			self.model.Xf=Xf
			self.model.Pf=Pf
			self.model.xi_est=xi_est
			self.model.kernel_weights_est=kernel_weights_est
			self.model.x0=Xb[0]
			self.model.P0=Pb[0]
			print it_count, " Kernel current estimate: ", self.model.kernel.weights, "xi", self.model.xi
			if it_count == max_it:
				keep_going = 0
			it_count += 1
		print "Elapsed time in seconds is", time.time()-t0

def dots(*args):
	lastItem = 1.
	for arg in args:
		lastItem = pb.dot(lastItem, arg)
	return lastItem
