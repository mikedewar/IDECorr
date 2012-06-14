from __future__ import division
import pylab as pb

class ActivationFunction():

	"""class defining the sigmoidal activation function .

	Arguments
	----------
	v0: float
		firing threshold, mV                  (Wendling, 2002, 6 mV), (Schiff, 2007, threshold = 0.24 (Heaviside))
	fmax: float
		maximum firing rate, spikes/s         (Wendling, 2002, 2*e0 = 5, or nu = 5
	varsigma: float
		slope of sigmoid, spikes/mV           (Wendling, 2002, 0.56 mV^-1)

	v: float
		Presynaptic potential

	Returns
	----------
	average firing rate
	"""	
	def __init__(self,v0,fmax,varsigma):
		
		self.v0 = v0
		self.fmax = fmax
		self.varsigma = varsigma
		
	def __call__(self,v):

		return float(self.fmax/(1.+pb.exp(self.varsigma*(self.v0-v))))


	def plot(self,plot_range):
		u=pb.linspace(-plot_range,plot_range,1000)
		z=pb.zeros_like(u)
		for i,j in enumerate(u):
			z[i]=self(j)
		pb.plot(u,z)
		pb.show()


