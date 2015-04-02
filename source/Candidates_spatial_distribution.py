#!/usr/bin/python


import numpy as np
import scipy as sc
import pylab as plt

GC       = sc.genfromtxt('Globular_Cluster.dat')
#data     = sc.genfromtxt('Candidata_Potenciales_V2_RAVE.dat')#'candidates_V2.dat')
data = sc.genfromtxt('Potenciales2_Candidates.dat')
data2 = sc.genfromtxt('Potenciales22_Candidates.dat')
fulldata = sc.genfromtxt('RAVE_Kinematics_D83.out')
SM       = sc.genfromtxt('Literature_Majewski_2012.dat')
GD       = sc.genfromtxt('Literature_DaCosta2008.dat')
WB       = sc.genfromtxt('Literature_WyliedeBoer2010.dat')
NS       = sc.genfromtxt('Literature_Schuster_2010.dat')

factor   = np.pi/180.

plt.subplot(2,1,1)
Nstar = 8.

plt.plot(NS[:,1],NS[:,2],'o',color='black')
plt.plot(SM[:,1],SM[:,2],'+',color='blue',ms=10,mew=2.)
plt.plot(GD[:,1],GD[:,2],'*',color='cyan',mec='cyan',ms=10,mew=2.)
plt.plot(WB[:,1],WB[:,2],'^',color='orange',mec='orange',ms=10,mew=2.,fillstyle='none')
plt.plot(data[:,1], data[:,2], 's',ms=Nstar,color='red',mew=2.,fillstyle='none')
plt.plot(data2[:,1], data2[:,2], 'D',ms=8,mew=2,color='red')
plt.plot(201.69683,-47.47958,'x',ms=17,mew=3.,color='green',fillstyle='none')
#plt.plot(GC[:,3],GC[:,4],'o',fillstyle='none',color='black')


# Kunder et al. (2014)
#plt.plot((260.,296.,296.,260.,260.),(-46.,-46.,-10.,-10.,-46.),'-',lw = 2., color ='gray')
plt.plot((70.,90.,90.,70.,70.),(-50.,-50.,-30.,-30.,-50.),'-',lw = 2., color ='gray')
plt.plot((152.,157.,157.,152.,152.),(-49.,-49.,-44.,-44.,-49.),'-',lw = 2., color ='gray')

plt.xlabel(r'$\alpha [^\circ]$',fontsize=20)
plt.ylabel(r'$\delta [^\circ]$',fontsize=20)
plt.tick_params(labelsize = 20)


plt.subplot(2,1,2)

maskl = (fulldata[:,32] < 180.)
plt.plot(fulldata[maskl,32],fulldata[maskl,30],'.',color='gray',alpha=0.02)
maskl = (fulldata[:,32] >= 180.)
plt.plot(fulldata[maskl,32]-360.,fulldata[maskl,30],'.',color='gray',alpha=0.02)

# Candidates RAVE

maskl = (data[:,18] < 180.)
plt.plot(data[maskl,18],data[maskl,140],'s',ms=Nstar,mew=2.,color='red',fillstyle='none')
maskl = (data[:,18] >= 180.)
plt.plot(data[maskl,18]-360.,data[maskl,140],'s',ms=Nstar,mew=2.,color='red',fillstyle='none')

maskl = (data2[:,18] < 180.)
plt.plot(data2[maskl,18],data2[maskl,140],'D',ms=8,mew=2,color='red')
maskl = (data2[:,18] >= 180.)
plt.plot(data2[maskl,18]-360.,data2[maskl,140],'D',ms=8,mew=2,color='red')

# Steven Majewski et al. (2012)

lSM , bSM        = SM[:,3]*factor, SM[:,4]*factor
Usun, Vsun, Wsun = 10.0, 5.23  , 7.17
Vclsr            = 220.
GRV              = SM[:,-1]
VrSM             = GRV - (Usun*np.cos(lSM)*np.cos(bSM) + (Vsun + Vclsr)*np.sin(lSM)*np.cos(bSM) + Wsun*np.sin(bSM))
Usun, Vsun, Wsun = np.abs(-11.1), np.abs(-12.24) , 7.25
Vclsr            = 239.0

GRV              = VrSM + Usun*np.cos(lSM)*np.cos(bSM) + (Vsun + Vclsr)*np.sin(lSM)*np.cos(bSM) + Wsun*np.sin(bSM)
Vb_star          = GRV / np.cos(bSM)
maskl = (SM[:,3] < 180.)
plt.plot(SM[maskl,3],Vb_star[maskl],'+',color='blue',ms=10,mew=2.)
maskl = (SM[:,3] >= 180.)
plt.plot(SM[maskl,3]-360.,Vb_star[maskl],'+',color='blue',ms=10,mew=2.)

# Da Costa et al. (2008)

lGD, bGD, VrGD = GD[:,3]*factor, GD[:,4]*factor, GD[:,5]

GRV              = VrGD + Usun*np.cos(lGD)*np.cos(bGD) + (Vsun + Vclsr)*np.sin(lGD)*np.cos(bGD) + Wsun*np.sin(bGD)
Vb_star          = GRV / np.cos(bGD)

plt.plot(GD[:,3]-360.,Vb_star,'*',color='cyan',mec='cyan',ms=10,mew=2.)

# Wylie-de Boer et al. (2010)


lWB, bWB, VrWB   = WB[:,-2]*factor, WB[:,-1]*factor, WB[:,5]
GRV              = VrWB + Usun*np.cos(lWB)*np.cos(bWB) + (Vsun + Vclsr)*np.sin(lWB)*np.cos(bWB) + Wsun*np.sin(bWB)
Vb_star          = GRV / np.cos(bWB)

maskl = (WB[:,-2] < 180.)
plt.plot(WB[maskl,-2],Vb_star[maskl],'^',color='orange',mec='orange',ms=10,mew=2.,fillstyle='none')
maskl = (WB[:,-2] >= 180.)
plt.plot(WB[maskl,-2]-360.,Vb_star[maskl],'^',color='orange',mec='orange',ms=10,mew=2.,fillstyle='none')

# Omega Centauri

l                = 309.10*factor
b                = 14.97*factor
Vr               = 232.1
Vclsr            = 239.0
Usun, Vsun, Wsun = np.abs(-11.1), np.abs(-12.24) , 7.25

GRV              = Vr + Usun*np.cos(l)*np.cos(b) + (Vsun + Vclsr)*np.sin(l)*np.cos(b) + Wsun*np.sin(b)
Vb_star          = GRV / np.cos(b)
plt.plot(309.10-360.,Vb_star,'x',ms=15,mew=3.,color='green',fillstyle='none')



plt.ylim(-500.,500.)
plt.xlim(180.,-180.)
plt.ylabel(r'$V_{GSR}/cos(b)$',fontsize=20.)
plt.xlabel('Galactic Longitude',fontsize=20.)
plt.tick_params(labelsize = 20)


plt.show()
