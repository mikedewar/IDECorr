from __future__ import division
import pylab as pb
from scipy import io
from scipy import signal
import os

#Read data
#~~~~~~~~~
#Y=io.loadmat('Y_50sec.mat')
Y=io.loadmat('data.mat')
Y=Y['Y'] 
#Transform observations in matrix forms
Data=pb.squeeze(Y)
Data.shape=Data.shape[0],pb.sqrt(Data.shape[1]),pb.sqrt(Data.shape[1])




AutoCorr=0
AutoCorrPlus=0
XCorr=0
NXCorr=0
for i in range(Data.shape[0]):
	AutoCorr+=signal.correlate2d(Data[i],Data[i],mode='full')


for i in range(Data.shape[0]-1):
	AutoCorrPlus+=signal.correlate2d(Data[i+1],Data[i+1],mode='full')
	XCorr+=signal.correlate2d(Data[i+1],Data[i],mode='full')
	NXCorr+=signal.correlate2d(Data[i],Data[i+1],mode='full')


#Calculate the mean
#~~~~~~~~~~~~~~~~~~

AutoCorr=AutoCorr/Data.shape[0]
AutoCorrPlus=AutoCorrPlus/(Data.shape[0]-1)
XCorr=XCorr/(Data.shape[0]-1)
NXCorr=NXCorr/(Data.shape[0]-1)

#save data
#~~~~~~~~~~
AutoCorr_temp={}
AutoCorrPlus_temp={}
XCorr_temp={}
NXCorr_temp={}

AutoCorr_temp['AutoCorr']=AutoCorr
AutoCorrPlus_temp['AutoCorrPlus']=AutoCorrPlus
XCorr_temp['XCorr']=XCorr
NXCorr_temp['NXCorr']=NXCorr

io.savemat('AutoCorr',AutoCorr_temp)
io.savemat('AutoCorrPlus',AutoCorrPlus_temp)
io.savemat('XCorr',XCorr_temp)
io.savemat('NXCorr',NXCorr_temp)






