from __future__ import division
import pylab as pb
from scipy import io
from scipy import signal
import os
import copy
from Bases2D import *
#Define parameters
#~~~~~~~~~~~~~~~~~~
#sampling properties

Fs = 1e3   #sampling rate                                       
Ts = 1/Fs   #sampling period, second

#Observation noise
ny=196
delta_s=1.5
varepsilon=0
dimension=2




#~~~~~~~~~~~~
#Reading Data
#~~~~~~~~~~~~~
AutoCorr_temp=io.loadmat('AutoCorr')
XCorr_temp=io.loadmat('XCorr')
AutoCorr=AutoCorr_temp['AutoCorr']
XCorr=XCorr_temp['XCorr']

#Divide correlation by number of observations
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
AutoCorr=AutoCorr/ny
XCorr=XCorr/ny

#Calculate Fourier transform of w(.)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#shift correlation to have the maximum at the top left corner( index 0,0)
ShiftedAutoCorr=pb.flipud(pb.fliplr(pb.fftshift(AutoCorr)))
ShiftedXCorr=pb.flipud(pb.fliplr(pb.fftshift(XCorr)))

#deep copy AutoCorrShift
ShiftedAutoCorrCopy=copy.deepcopy(ShiftedAutoCorr)
#subtract observation noise variance from the maximum,(from index 0,0)
ShiftedAutoCorrCopy[0,0]=ShiftedAutoCorrCopy[0,0]-varepsilon
Denum=pb.fft2(ShiftedAutoCorrCopy)

#calculate Fourier Transform of the kernel
FFTShiftedXCorr=pb.fft2(ShiftedXCorr)
W=pb.conjugate((FFTShiftedXCorr/Denum))
#W=(FFTShiftedXCorr/Denum)

#find the kernel in spatial domain
w=pb.real(pb.ifft2(W))/(delta_s**2)
w=pb.ifftshift(w)
w=pb.roll(w,-1,axis=0)
w=pb.roll(w,-1,axis=1)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Define the true kernel

#centers
#psi0_center=pb.matrix([[-0.5],[-0.5]]);psi1_center=pb.matrix([[0],[0]]);psi2_center=pb.matrix([[0.5],[0.5]])

psi0_center=pb.matrix([[0],[0]]);psi1_center=pb.matrix([[0],[0]]);psi2_center=pb.matrix([[0],[0]])
#widths
psi0_width=1.8**2;psi1_width=2.4**2;psi2_width=6.**2
#psi0_width=2.4**2;psi1_width=2.4**2;psi2_width=2.4**2

#weights
#psi0_weight=-.07;psi1_weight=.07;psi2_weight=.02
psi0_weight=.5;psi1_weight=-.3;psi2_weight=0

#define each kernel basis functions
psi0=basis(psi0_center,psi0_width,dimension)
psi1=basis(psi1_center,psi1_width,dimension)
psi2=basis(psi2_center,psi2_width,dimension)
#Define space
#~~~~~~~~~~~~~
s=pb.linspace(-20,20,len(w.diagonal()))
#s=pb.linspace(-20,20,500)
s1,s2=pb.meshgrid(s,s)
Psi=psi0_weight*psi0(s1,s2)+psi1_weight*psi1(s1,s2)+psi2_weight*psi2(s1,s2)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#plot the results
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ext=20
#Mexican hat
#~~~~~~~~~~~
vmin=-0.03;vmax=0.2
fig1=pb.figure(figsize=(4,2.39))
ax11=fig1.add_subplot(121)
cax11=ax11.imshow(Psi,vmin=vmin,vmax=vmax,extent=[-ext,ext,-ext,ext],interpolation='nearest')#,cmap=pb.cm.hot
cbar11=fig1.colorbar(cax11,orientation='horizontal',shrink=0.6,ticks=[-0.03,.2])#
cbar11.ax.xaxis.set_ticks_position('top')
ax11.set_xticks([-ext,0,ext])
ax11.set_yticks([0,ext])
ax11.set_xlabel('Space',fontname='Arial',fontsize=10)
ax11.set_ylabel('Space',fontname='Arial',fontsize=10)

#ax11.text(-15,-2.9,'Space',rotation=90,fontname='Arial',fontsize=10)
#ax11.text(-3,-15,'Space',fontname='Arial',fontsize=10)
#ax11.text(12,-2.9,'Space',rotation=90,fontname='Arial',fontsize=10)
#ax11.text(22.9,-15,'Space',fontname='Arial',fontsize=10)
ax12=fig1.add_subplot(122)
cax12=ax12.imshow(w,vmin=vmin,vmax=vmax,extent=[-ext,ext,-ext,ext],interpolation='nearest')#,cmap=pb.cm.hot
cbar12=fig1.colorbar(cax12,orientation='horizontal',shrink=.6,ticks=[-0.03,.2])
cbar12.ax.xaxis.set_ticks_position('top')
ax12.set_xlabel('Space',fontname='Arial',fontsize=10)
ax12.set_xticks([-ext,0,ext])
ax12.set_yticks([0,ext])

#psitioning the colorbars
sub1=fig1.get_axes()[0]
sub1_cbar=fig1.get_axes()[1]
sub2=fig1.get_axes()[-2]
sub2_cbar=fig1.get_axes()[-1]

sub1_cbar.set_position([0.17, 0.23, 0.22, 0.65])
sub2_cbar.set_position([0.65, 0.23, 0.22, 0.65])
fig1.subplots_adjust(left=0.07, bottom=0.23, right=0.97, top=0.79,wspace=0.01, hspace=0.56)
fig1.savefig('KernelEstimation2d.pdf',format='pdf', dpi=300) 

fig2=pb.figure(figsize=(3.27,2))
fig2.subplots_adjust(left=0.24, bottom=0.23, right=0.94, top=0.94,wspace=0.2, hspace=0.2)
ax2=fig2.add_subplot(111)
ax2.plot(s,Psi.diagonal(),'k',label=r'$k$')
ax2.plot(s,w.diagonal(),'k:',label=r'$\hat{k}$')
leg=ax2.legend(frameon=False,prop={'size':'medium'})

ax2.set_xticks([-20,0,20])
ax2.set_yticks([-0.05,0.25])
ax2.set_xlabel('Space',fontsize=10,fontname='Arial')
ax2.text(-25,.14,'Amplitude',fontsize=10,fontname='Arial',rotation='vertical')
fig2.savefig('KernelSupportEstimation1d.pdf',format='pdf',bbox_inches='tight', dpi=300) 

pb.show()


