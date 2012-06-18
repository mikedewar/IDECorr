from __future__ import division
import pylab as pb
import numpy as np
import scipy as sp
from Bases2D import *
from scipy import signal
from scipy import io
# spatial parameters
# ~~~~~~~~~~~
Delta = 0.5;                          # space step for the spatial discretisation
Delta_squared = Delta**2
SpaceMaxPeriodicField = 30                   # maximum space in mm
SpaceMinPeriodicField = -30         # minimum space in mm
NPointsInPeriodicField = int((SpaceMaxPeriodicField-SpaceMinPeriodicField)/Delta+1)
NPointsInField = int((NPointsInPeriodicField-1)/3 + 1)
r = pb.linspace(SpaceMinPeriodicField/3,SpaceMaxPeriodicField/3,NPointsInField)      # define space
r_periodic=pb.linspace(SpaceMinPeriodicField,SpaceMaxPeriodicField,NPointsInPeriodicField) #define periodic space
r1,r2=pb.meshgrid(r_periodic,r_periodic)  


# disturbance paramters
# ~~~~~~~~~~~~~
gamma_width=1.3**2        # parameter for covariance of disturbance
gamma_weight = 1            # variance of disturbance
dimension=2
# ~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~
FirstThirdEnd = int((NPointsInPeriodicField-1)/3)-1 #we leave out the center, FirstThird ends at 39 [0,39] equals to 40 points index wise [0:40]
SecondThirdEnd = 2*FirstThirdEnd+2 #second third is [-10,10] equals to 41 to 80 points, index wise: 41:81

mm=0
Sigma_gamma = pb.zeros((NPointsInField**2,NPointsInField**2))   #create disturbance covariance matrix
for n in range(0,NPointsInField):
    for nn in range(0,NPointsInField):

	gamma=basis(pb.matrix([r[n],r[nn]]),gamma_width,dimension)

        temp = gamma_weight*gamma(r1,r2)

	bottomleft_temp=pb.hstack((pb.zeros((FirstThirdEnd+1,1)),temp[0:FirstThirdEnd+1,0:FirstThirdEnd+1]))
	bottomleft=pb.vstack((pb.zeros((1,FirstThirdEnd+2)),bottomleft_temp))

	bottom=pb.vstack((pb.zeros((1,FirstThirdEnd+2)),temp[0:FirstThirdEnd+1,FirstThirdEnd+1:SecondThirdEnd+1]))

	bottomright_temp=pb.hstack((temp[0:FirstThirdEnd+1,SecondThirdEnd+1:],pb.zeros((FirstThirdEnd+1,1))))
	bottomright=pb.vstack((pb.zeros((1,FirstThirdEnd+2)),bottomright_temp))

	left=pb.hstack((pb.zeros((FirstThirdEnd+2,1)),temp[FirstThirdEnd+1:SecondThirdEnd+1,0:FirstThirdEnd+1]))
	middle=temp[FirstThirdEnd+1:SecondThirdEnd+1,FirstThirdEnd+1:SecondThirdEnd+1]
	right=pb.hstack((temp[FirstThirdEnd+1:SecondThirdEnd+1,SecondThirdEnd+1:],pb.zeros((FirstThirdEnd+2,1))))


	topleft_temp=pb.vstack((temp[SecondThirdEnd+1:,0:FirstThirdEnd+1],pb.zeros((1,FirstThirdEnd+1))))
	topleft=pb.hstack((pb.zeros((FirstThirdEnd+2,1)),topleft_temp))

	top=pb.vstack((temp[SecondThirdEnd+1:,FirstThirdEnd+1:SecondThirdEnd+1],pb.zeros((1,FirstThirdEnd+2))))


	topright_temp=pb.hstack((temp[SecondThirdEnd+1:,SecondThirdEnd+1:],pb.zeros((FirstThirdEnd+1,1))))
	topright=pb.vstack((topright_temp,pb.zeros((1,FirstThirdEnd+2))))



	temp2=middle+topleft+top+topright+left+right+bottom+bottomleft+bottomright
	
	
	Sigma_gamma[:,mm] = temp2.ravel('F');
        mm=mm+1;

data=dict(Sigma_gamma=Sigma_gamma)
io.savemat('periodic_sigma_gamma',data)




