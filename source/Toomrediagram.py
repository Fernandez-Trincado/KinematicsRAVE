#!/usr/bin/python


import numpy as np
import scipy as sc
import pylab as plt

sim  = sc.genfromtxt('Simu_m1407_selHRV-9-12.dat1_D83.out')
#rave = sc.genfromtxt('Candidata_Potenciales_V2_RAVE.dat')
rave  = sc.genfromtxt('Potenciales2_Candidates.dat')
rave2 = sc.genfromtxt('Potenciales22_Candidates.dat')
fullra = sc.genfromtxt('RAVE_Kinematics_D83.out')
add    = sc.genfromtxt('RAVE_Kinematics_D83_cuts.out')
#sch    = sc.genfromtxt('Literature_Schuster_2010.dat')

mask_  = (sim[:,22]>=240. ) & (sim[:,22]<=360.) & (sim[:,23]>=-60.) & (sim[:,23]<=60.)
sim    = sim[mask_] 

mask_2 = (fullra[:,32]>=240. ) & (fullra[:,32]<=360.) & (fullra[:,33]>=-60.) & (fullra[:,33]<=60.)
fullra = fullra[mask_2]


# Simulations

thin     = (sim[:,16] <= 7.)
Ythick   = (sim[:,16] == 8.)
#Othick   = (sim[:,16] == 11.)
halo     = (sim[:,16] == 9.)
bar      = (sim[:,16] == 10.)

UUsim    = sim[:,8] 
VVsim    = sim[:,9]
WWsim    = sim[:,10]

plt.subplot(1,1,1)

Toomresim = np.sqrt((UUsim*UUsim) + (WWsim*WWsim))

plt.plot(VVsim[Ythick], Toomresim[Ythick], '.' ,fillstyle='none' ,color='orange')
plt.plot(VVsim[thin]  , Toomresim[thin]  , '.' ,color='green')
#plt.plot(VVsim[Othick], Toomresim[Othick], '.' ,color='magenta')
plt.plot(VVsim[halo]  , Toomresim[halo]  , '.' ,color='cyan')

#UUadd = add[:,25]
#VVadd = add[:,26] + 239.
#WWadd = add[:,27]
#
#Toomreadd = np.sqrt((UUadd*UUadd) + (WWadd*WWadd))
#plt.plot(VVadd, Toomreadd, '*',mec='green',color='green')

# Candidates RAVE

UUrave    = rave[:,135]
VVrave    = rave[:,142] #+ 239.
WWrave    = rave[:,137]

Toomrerave = np.sqrt((UUrave*UUrave) + (WWrave*WWrave))
plt.plot(VVrave, Toomrerave, 's',ms=8,mew=2, fillstyle='none' ,color='red')
plt.xlim(-600.,200.)
plt.ylim(0.,700.)
plt.ylabel(r'$(U_{LSR}^2 + W_{LSR}^2)^{1/2}$',fontsize=20)
plt.xlabel(r'$V_{LSR}$ (km/s)',fontsize=20)
plt.axvline(x=-239.,ls='--',lw=2.,color='black')
plt.text(-550.,550.,'Retrograde',color='black',fontsize=20)
plt.text(-100.,550.,'Prograde',color='black',fontsize=20)


UUrave    = rave2[:,135]
VVrave    = rave2[:,142] #+ 239.
WWrave    = rave2[:,137]

Toomrerave = np.sqrt((UUrave*UUrave) + (WWrave*WWrave))
plt.plot(VVrave, Toomrerave, 'D',ms=8,mew=2,color='red')

# Schuster 

#UUsch    = sch[:,-3]
#VVsch    = sch[:,-2]
#WWsch    = sch[:,-1]
#
#Toomresch = np.sqrt((UUsch*UUsch) + (WWsch*WWsch)) 
#plt.plot(VVsch, Toomresch, '*',fillstyle='none',ms=12,mew=1, color='black')


###################################################################################
#plt.subplot(1,2,2)
#
#UUadd = add[:,25]
#VVadd = add[:,26] + 239.
#WWadd = add[:,27]
#
#Toomreadd = np.sqrt((UUadd*UUadd) + (WWadd*WWadd))
#plt.plot(VVadd, Toomreadd, '.',mec='gray',color='gray')
#
##UUf = fullra[:,25]
##VVf = fullra[:,26] + 239.
##WWf = fullra[:,27]
##
##Toomref = np.sqrt((UUf*UUf) + (WWf*WWf))
##plt.plot(VVf, Toomref, '.',alpha=0.06,color='gray')
#plt.xlim(-600.,200.)
#plt.ylim(0.,700.)
#plt.ylabel(r'$(U_{LSR}^2 + W_{LSR}^2)^{1/2}$ (km/s)',fontsize=20)
#plt.xlabel(r'$V_{LSR}$ (km/s)',fontsize=20)
#plt.text(-450.,550.,'RAVE survey',color='black',fontsize=20)
#

#plt.subplot(1,3,3)
#
#Othick   = (sim[:,16] != 11.)
#
#plt.plot(VVsim[Othick], Toomresim[Othick], '.',alpha=0.06,color='gray')
#plt.xlim(-600.,200.)
#plt.ylim(0.,700.)
#plt.ylabel(r'$(U_{LSR}^2 + W_{LSR}^2)^{1/2}$ (km/s)',fontsize=20)
#plt.xlabel(r'$V_{LSR}$ (km/s)',fontsize=20)
#plt.text(-450.,550.,'Besancon Galaxy Model', color='black',fontsize=20)



plt.show()
