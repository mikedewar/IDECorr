'''This is to generate IDE components: kernal and the field'''


class Kernel():

	"""class defining the Connectivity kernel of the brain.

	Arguments
	----------

	Psi:list
		list of connectivity kernel basis functions; must be of class Bases
	weights: list
		list of connectivity kernel weights	


	Attributes
	----------
	evaluate: 
		evaluate the kernel at a given spatial location
	"""
	
	def __init__(self,Psi,weights):
		
		self.Psi=Psi
		self.weights = weights



class Field():	

	"""class defining the field.

	Arguments
	----------

	Psi:list
		list of field basis functions; must be of class Bases

	Attributes
	----------
	evaluate: 
		evaluate the kernel at a given spatial location
	"""

	def __init__(self, Phi):

		self.Phi=Phi


