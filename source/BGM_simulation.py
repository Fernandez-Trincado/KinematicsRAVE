#!/usr/bin/python


import numpy as np
import scipy as sc
import pylab as plt

data = sc.genfromtxt('Simu_m1407_selHRV-9-12.dat1_D83.out')
 
l    = data[:,22]
b    = data[:,23]
Vr   = data[:,7]
GRV  = data[:,44]
Vb   = data[:,43]

thin     = (data[:,16] <= 7.)
Ythick   = (data[:,16] == 8.)
#Othick   = (data[:,16] == 11.)
halo     = (data[:,16] == 9.)
#bar      = (data[:,16] == 10.)


# Velocities
mask = (l>=240.) & (l<= 360.) & (b>=-60.) & (b <= 60.) & (Vr >= 95. ) & (Vr <= 300. ) 
thin   = thin[mask]
Ythick = Ythick[mask]
halo   = halo[mask]

Common = GRV[mask]

#plt.hist(Common[thin]  ,lw=2.,histtype='stepfilled',color='orange',edgecolor="orange",fill=False,range=(-500.,500.),bins=1000./50.)
plt.hist(Common[Ythick],lw=2.,histtype='stepfilled',color='gray',edgecolor="gray",fill=False,range=(-500.,500.),bins=1000./50.)
plt.hist(Common[halo]  ,lw=2.,histtype='stepfilled',color='blue',edgecolor="blue",fill=False,range=(-500.,500.),bins=1000./50.)


plt.show()


