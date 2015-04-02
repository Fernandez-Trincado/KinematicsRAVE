#!/usr/bin/python


import numpy as np
import pylab as plt
import scipy as sc


data = sc.genfromtxt('Simu_m1407_selHRV-9-12.dat1_D83.out')

Vr   = data[:,7]
l    = data[:,22]
b    = data[:,23]

mask_ = (l>=240.) & (l<= 360.) & (b>=-60.) & (b <= 60.) & (Vr >= 95. ) & (Vr <= 300. )
data  = data[mask_]

MH  = data[:,21]
RV  = data[:,7]
GRV = data[:,44]
VbS = GRV#data[:,-6]+239.#data[:,43]  

thin     = (data[:,16] <= 7.)
Ythick   = (data[:,16] == 8.)
#Othick   = (data[:,16] == 11.)
halo     = (data[:,16] == 9.)
bar      = (data[:,16] == 10.)

plt.subplot(1,2,1)

plt.plot(VbS[thin]   , MH[thin]   , '.' ,color='gray')                     # Thin disc
plt.plot(VbS[Ythick] , MH[Ythick] , '.' ,fillstyle='none',color='orange')   # Ythick disc
#plt.plot(RV[Othick] , MH[Othick] , 'x' ,fillstyle='none',color='orange')  # Ythick disc
plt.plot(VbS[halo]   , MH[halo]   , '.' ,mec='cyan',color='cyan')    # halo
plt.plot(VbS[bar]    , MH[bar]    , '.' ,fillstyle='none',color='magenta') # bar
plt.xlim(-400.,400.)
plt.ylim(-3.0,1.0)

# Data RAVE

#RAVE = sc.genfromtxt('Candidata_Potenciales_V2_RAVE.dat')
RAVE = sc.genfromtxt('Potenciales2_Candidates.dat')
RAVE2 = sc.genfromtxt('Potenciales22_Candidates.dat')


alpha  = RAVE[:,94]
FeH    = RAVE[:,83]
falpha = 10**(alpha)
old_mH = RAVE[:,60]
RVRAVE = RAVE[:,5]
GRVR   = RAVE[:,141]
VbR    = GRVR#RAVE[:,-7]#RAVE[:,131]

MH_RAVE = FeH + (np.log( 0.638*falpha + 0.362 )/np.log(10))

#print MH_RAVE, old_mH

plt.plot(VbR, MH_RAVE, 's', ms=8,mew=2, fillstyle='none', color='red')


alpha  = RAVE2[:,94]
FeH    = RAVE2[:,83]
falpha = 10**(alpha)
old_mH = RAVE2[:,60]
RVRAVE = RAVE2[:,5]
GRVR   = RAVE2[:,141]
VbR    = GRVR#RAVE[:,-7]#RAVE[:,131]

MH_RAVE = FeH + (np.log( 0.638*falpha + 0.362 )/np.log(10))

#print MH_RAVE, old_mH

plt.plot(VbR, MH_RAVE, 'D',ms=8,mew=2,color='red')


#RAVE = sc.genfromtxt('Potenciales3_Candidates.dat')
#
#alpha  = RAVE[:,94]
#FeH    = RAVE[:,83]
#falpha = 10**(alpha)
#old_mH = RAVE[:,60]
#RVRAVE = RAVE[:,5]
#GRVR   = RAVE[:,141]
#VbR    = GRVR#RAVE[:,-7]#RAVE[:,131]
#
#MH_RAVE = FeH + (np.log( 0.638*falpha + 0.362 )/np.log(10))
#
##print MH_RAVE, old_mH
#
#plt.plot(VbR, MH_RAVE, 's', ms=8,mew=2, fillstyle='none', color='red')


#plt.axvline(x=220.,ls='--',lw=2,color='black')
#plt.axvline(x=280.,ls='--',lw=2,color='black')
plt.tick_params(labelsize = 20)
plt.xlabel(r'$V_{GSR}$ (km/s)',fontsize=20)
#plt.xlabel(r'$V_{GRS}/cos(b)$',fontsize=20)
plt.ylabel('[M/H] dex',fontsize=20)






data = sc.genfromtxt('Simu_m1407_selHRV-9-12.dat1_D83.out')

Vr   = data[:,7]
l    = data[:,22]
b    = data[:,23]

mask_ = (l>=240.) & (l<= 360.) & (b>=-60.) & (b <= 60.) #& (Vr >= 95. ) & (Vr <= 300. )
data  = data[mask_]

MH  = data[:,21]
RV  = data[:,7]
GRV = data[:,44]
VbS = GRV#data[:,-6]+239.#data[:,43]  

thin     = (data[:,16] <= 7.)
Ythick   = (data[:,16] == 8.)
#Othick   = (data[:,16] == 11.)
halo     = (data[:,16] == 9.)
bar      = (data[:,16] == 10.)

plt.subplot(1,2,2)

plt.plot(VbS[thin]   , MH[thin]   , '.' ,color='gray')                     # Thin disc
plt.plot(VbS[Ythick] , MH[Ythick] , '.' ,fillstyle='none',color='orange')   # Ythick disc
#plt.plot(RV[Othick] , MH[Othick] , 'x' ,fillstyle='none',color='orange')  # Ythick disc
plt.plot(VbS[halo]   , MH[halo]   , '.' ,mec='cyan',color='cyan')    # halo
plt.plot(VbS[bar]    , MH[bar]    , '.' ,fillstyle='none',color='magenta') # bar

# Data RAVE

#RAVE = sc.genfromtxt('Candidata_Potenciales_V2_RAVE.dat')
RAVE = sc.genfromtxt('Potenciales2_Candidates.dat')
RAVE2 = sc.genfromtxt('Potenciales22_Candidates.dat')


alpha  = RAVE[:,94]
FeH    = RAVE[:,83]
falpha = 10**(alpha)
old_mH = RAVE[:,60]
RVRAVE = RAVE[:,5]
GRVR   = RAVE[:,141]
VbR    = GRVR# RAVE[:,-7]#RAVE[:,131]

MH_RAVE = FeH + (np.log( 0.638*falpha + 0.362 )/np.log(10))

#print MH_RAVE, old_mH

plt.plot(VbR, MH_RAVE, 's', ms=8,mew=2, fillstyle='none', color='red')

alpha  = RAVE2[:,94]
FeH    = RAVE2[:,83]
falpha = 10**(alpha)
old_mH = RAVE2[:,60]
RVRAVE = RAVE2[:,5]
GRVR   = RAVE2[:,141]
VbR    = GRVR# RAVE[:,-7]#RAVE[:,131]

MH_RAVE = FeH + (np.log( 0.638*falpha + 0.362 )/np.log(10))

#print MH_RAVE, old_mH

plt.plot(VbR, MH_RAVE, 'D',ms=8,mew=2,color='red')



#RAVE = sc.genfromtxt('Potenciales3_Candidates.dat')
#
#alpha  = RAVE[:,94]
#FeH    = RAVE[:,83]
#falpha = 10**(alpha)
#old_mH = RAVE[:,60]
#RVRAVE = RAVE[:,5]
#GRVR   = RAVE[:,141]
#VbR    = GRVR#RAVE[:,-7]#RAVE[:,131]
#
#MH_RAVE = FeH + (np.log( 0.638*falpha + 0.362 )/np.log(10))
#
##print MH_RAVE, old_mH
#
#plt.plot(VbR, MH_RAVE, 's', ms=8,mew=2, fillstyle='none', color='red')

plt.plot(VbR, MH_RAVE, 'D',ms=8,mew=2,color='red')


#plt.axvline(x=220.,ls='--',lw=2,color='black')
#plt.axvline(x=280.,ls='--',lw=2,color='black')
plt.tick_params(labelsize = 20)
plt.xlabel(r'$V_{GSR}$ (km/s)',fontsize=20)
#plt.xlabel(r'$V_{GRS}/cos(b)$',fontsize=20)
plt.ylabel('[M/H] dex',fontsize=20)
plt.xlim(-400.,400.)
plt.ylim(-3.0,1.0)




plt.show()


