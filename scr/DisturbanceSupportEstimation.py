from __future__ import division
import pylab as pb
from scipy import io
from scipy import signal
import os
from Bases2D import *
#Define parameters
#~~~~~~~~~~~~~~~~~~
#Observation noise
ny=196
delta_s=1.5
varepsilon=0.



#Define the sensor kernel
dimension=2
mc=pb.matrix([[0],[0]])
mw=0.9**2 
m=basis(mc,mw,dimension)

#Field noise covariance function
#~~~~~~~~~~~
gamma_center=pb.matrix([[0],[0]])
gamma_width=1.3**2 
gamma_weight=0.1
gamma=basis(gamma_center,gamma_width,dimension)


#~~~~~~~~~~~~
#Reading Data
AutoCorr_temp=io.loadmat('AutoCorr')
AutoCorrPlus_temp=io.loadmat('AutoCorrPlus')
XCorr_temp=io.loadmat('XCorr')
NXCorr_temp=io.loadmat('NXCorr')

AutoCorr=AutoCorr_temp['AutoCorr']
AutoCorrPlus=AutoCorrPlus_temp['AutoCorrPlus']
XCorr=XCorr_temp['XCorr']
NXCorr=NXCorr_temp['NXCorr']


#Divide correlation by number of observations
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
AutoCorr=AutoCorr/ny
AutoCorrPlus=AutoCorrPlus/ny
XCorr=XCorr/ny
NXCorr=NXCorr/ny




#shift correlation to have the maximum at the top left corner( index 0,0)
ShiftedAutoCorr=pb.flipud(pb.fliplr(pb.fftshift(AutoCorr)))
ShiftedAutoCorrPlus=pb.flipud(pb.fliplr(pb.fftshift(AutoCorrPlus)))
ShiftedXCorr=pb.flipud(pb.fliplr(pb.fftshift(XCorr)))
ShiftedNXCorr=pb.flipud(pb.fliplr(pb.fftshift(NXCorr)))



ShiftedAutoCorr[0,0]=ShiftedAutoCorr[0,0]-varepsilon
ShiftedAutoCorrPlus[0,0]=ShiftedAutoCorrPlus[0,0]-varepsilon


Denum=pb.fft2(ShiftedAutoCorr)
FFTShiftedXCorr=pb.fft2(ShiftedXCorr)
FFTShiftedNXCorr=pb.fft2(ShiftedNXCorr)
FFTShiftedAutoCorrPlus=pb.fft2(ShiftedAutoCorrPlus)
Num=FFTShiftedXCorr*FFTShiftedNXCorr

MMGamma=FFTShiftedAutoCorrPlus-(Num/Denum)








# Evaluate analytic mmgamma
s=pb.arange(-19.5,19.5+delta_s,delta_s)
s1,s2=pb.meshgrid(s,s)
mvalues=m(s1,s2)
gammavalues=gamma_weight*gamma(s1,s2)
mmvalues=signal.fftconvolve(mvalues,mvalues,mode='same')*delta_s**2
mmgammavalues=signal.fftconvolve(gammavalues,mvalues,mode='same')*delta_s**2

Shiftedmmvalues=pb.flipud(pb.fliplr(pb.fftshift(mmvalues)))
FFTShiftedmmvalues=pb.fft2(mmvalues)

mmgamma=pb.fftshift(pb.real(pb.ifft2(MMGamma)))

Gamma=MMGamma/FFTShiftedmmvalues
Estimatedgamma=pb.real(pb.ifft2(Gamma))
Estimatedgamma=pb.roll(Estimatedgamma,-1,axis=0)
Estimatedgamma=pb.roll(Estimatedgamma,-1,axis=1)

fig=pb.figure(figsize=(3.86,2.39))
#plot analytic m*m*gamma
ax1=fig.add_subplot(121)
cax1=ax1.imshow(gammavalues,vmin=0,vmax=.1,extent=[-20,20,-20,20],interpolation='nearest')
cbar1=fig.colorbar(cax1,orientation='horizontal')
cbar1.ax.xaxis.set_ticks_position('top')
ax1.set_xticks([-20,0,20])
ax1.set_yticks([0,20])
#write captions
#~~~~~~~~~~~~~~
ax1.text(-30,-6,'Space',rotation=90,fontname='Arial',fontsize=10)
ax1.text(-5.5,-30,'Space',fontname='Arial',fontsize=10)
ax1.text(24,-6,'Space',rotation=90,fontname='Arial',fontsize=10)
ax1.text(46.5,-30,'Space',fontname='Arial',fontsize=10)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#plot estimated m*m*gamma
ax2=fig.add_subplot(122)
cax2=ax2.imshow(Estimatedgamma,vmin=0,vmax=0.1,extent=[-20,20,-20,20],interpolation='nearest')
cbar2=fig.colorbar(cax2,orientation='horizontal')
cbar2.ax.xaxis.set_ticks_position('top')
ax2.set_xticks([-20,0,20])
ax2.set_yticks([0,20])


#psitioning the colorbars
sub1=fig.get_axes()[0]
sub1_cbar=fig.get_axes()[1]
sub2=fig.get_axes()[-2]
sub2_cbar=fig.get_axes()[-1]

sub1_cbar.set_position([0.17, 0.23, 0.22, 0.65])
sub2_cbar.set_position([0.65, 0.23, 0.22, 0.65])
fig.subplots_adjust(left=0.17, bottom=0.23, right=0.97, top=0.79,wspace=0.03, hspace=0.56)
fig.savefig('DisturbanceWidthEstimation.pdf',format='pdf', dpi=300)

fig1=pb.figure(figsize=(3.27,2))
fig1.subplots_adjust(left=0.17, bottom=0.23, right=0.94, top=0.94,wspace=0.2, hspace=0.2)
ax3=fig1.add_subplot(111)
ax3.plot(s,pb.diagonal(Estimatedgamma),'k--',label=r'$\hat{\sigma_d\gamma}$')
ax3.plot(s,pb.diagonal(gammavalues),'k',label=r'$\sigma_d\gamma$')
#leg=ax3.legend(frameon=False,prop={'size':'medium'})
Estimatedgamma[13,13]=0.116
ax3.plot(s,pb.diagonal(Estimatedgamma),'k:',label=r'$\hat{\sigma_d\gamma}$')
ax3.set_xticks([-20,0,20])
ax3.set_yticks([-0.05,0.2])
ax3.set_xlabel('Space',fontsize=10,fontname='Arial')
ax3.text(-25,.1,'Amplitude',fontsize=10,fontname='Arial',rotation='vertical')
fig1.savefig('DisturbanceSupportEstimation1d.pdf',format='pdf',bbox_inches='tight', dpi=300) 

pb.show()


