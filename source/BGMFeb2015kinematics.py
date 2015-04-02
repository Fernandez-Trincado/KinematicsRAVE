#!/usr/bin/python
# Created by: J. G. Fernandez-Trincado
# Date 05/02/2015
# Last update: 28/02/2015i
# Computing U, V, W (km/s)

import numpy as np
import scipy as sc
import pylab as plt

file_  = sc.genfromtxt('Simu_m1407_selHRV-9-12.dat1',dtype=str) # Input file
factor = (np.pi/180.)
R_sun  = 8.0#8.3      # Brunthaler et al. (2011)
k      = 4.741

#output = open('Simu_m1407_selHRV-9-12.dat1.out','a')
#output.write('# RA_  DEC_  l b X  Y  Z  R_gal U V W VR  Vphi  Vb_star  GRV \n')

for i in np.arange(len(file_[:,0])):

    Vradial     = float(file_[i,7])
    D           = float(file_[i,26])
    l_,b_       = float(file_[i,22]), float(file_[i,23])
    RA_,DEC_    = float(file_[i,24]), float(file_[i,25])
    pmRA, pmDEC = float(file_[i,5]), float(file_[i,6])
    
    RAG, DECG   = 192.85948*factor, 27.12825*factor
    l, b        = l_*factor, b_*factor
    RA, DEC     = RA_*factor, DEC_*factor
    C1          = np.sin(DECG)*np.cos(DEC)-np.cos(DECG)*np.sin(DEC)*np.cos(RA-RAG)
    C2          = np.cos(DECG)*np.sin(RA-RAG)
    C           = np.matrix([[C1,C2],[-C2,C1]])
    pm          = np.matrix([[pmRA],[pmDEC]])
    pmG         = (1./np.cos(b))*C*pm
	

# Calculating Proper motion in galactic coordinates (mu_l, mu_b) from Radoslaw Polesky (2013)

    pml, pmb = np.cos(b)*pmG[0,0], pmG[1,0]

# Velocities (Vl, Vb) in galactic coordinates from Radoslaw Polesky (2013)

    Vl, Vb = D*k*np.cos(b)*pmG[0,0], D*k*pmG[1,0]

# Calculating X, Y and Z from Bond et al. (2010), ApJ, 716:1-29

    X     = R_sun-D*np.cos(l)*np.cos(b)
    Y     = -D*np.sin(l)*np.cos(b)
    Z     = D*np.sin(b)
    R_gal = np.sqrt((X*X)+(Y*Y)+(Z*Z)) # 3-dimension
    R     = np.sqrt((X*X)+(Y*Y))       # 2-dimension 

# Calculating VX, VY and VZ from Bond et al. (2010), ApJ, 716:1-29

    VX = -Vradial*np.cos(l)*np.cos(b) + Vb*np.cos(l)*np.sin(b) + Vl*np.sin(l)
    VY = -Vradial*np.sin(l)*np.cos(b) + Vb*np.sin(l)*np.sin(b) - Vl*np.cos(l)
    VZ =  Vradial*np.sin(b) + Vb*np.cos(b)

# Calculating U, V and W from Binney & Merrifield (1998)

#    UU, VV, WW = -VX, -VY , VZ # System Binney

# Calculating VR, Vphi ; Bond et al. (2010) 
# X axis is oriented toward l=180 degree, and the Y axis is oriented toward l = 270 degree, the disk rotates toward l = 90 degree.

    Vclsr            = 0 #239.0 # Brunthaler et al. (2011)   
    Ulsr, Vlsr, Wlsr =  0.,   Vclsr,   0.      
#	Usun, Vsun, Wsun = -10.0,  -5.23,   7.17   # Dehnen & Binney (1998) and Hogg et al. (2005)
#	Usun, Vsun, Wsun = -11.1, -12.24 , 7.25     # Brunthaler et al. (2011) and Schonrich et al. (2010)
    Usun, Vsun, Wsun = -11.0,    -12.0,     7.0

    VX   = VX + (Usun - Ulsr)
    VY   = VY + (Vsun - Vlsr)
    VZ   = VZ + (Wsun - Wlsr)

    UU, VV, WW = -VX, -VY , VZ

    print UU, VV, WW, file_[i,8], file_[i,9], file_[i,10]


# Our system    

#    UU, VV, WW = VX, VY , VZ
# Note that, in our adopted system, the disk has a prograde rotation Vphi = -220., retrograde rotation is indicated by Vphi > 0. Stars with VR > 0 move away from the Galactic center, and stars with Vz > 0 move toward the North Galactic Pole.	

    VR   = (VX*(X/R))+(VY*(Y/R))
    Vphi = -(VX*(Y/R))+(VY*(X/R))

# Calculating VLSR

    GRV       = Vradial + Usun*np.cos(l)*np.cos(b) + (Vsun + Vclsr)*np.sin(l)*np.cos(b) + Wsun*np.sin(b)
    Vb_star   = GRV / np.cos(b)
	
#	output.write(str(RA_)+' '+str(DEC_)+' '+' '+str(l)+' '+str(b)+' '+str(X)+' '+str(Y)+' '+str(Z)+' '+str(R_gal)+' '+str(UU)+' '+str(VV)+' '+str(WW)+' '+str(VR)+' '+str(Vphi)+' '+str(Vb_star)+' '+str(GRV)+' \n')
 
#output.close()
