#!/usr/bin/python
# Created by: J. G. Fernandez-Trincado
# Date 05/02/2015
# Last update: 28/02/2015i
# Computing U, V, W (km/s)

import numpy as np
import scipy as sc
import pylab as plt

#file_  = sc.genfromtxt('Simu_m1407_selHRV-9-12.dat1',dtype=str) # Input file
file_  = sc.genfromtxt('RAVE_DR4.dat',dtype=str) # Input file
factor = (np.pi/180.)
k      = 4.741

#R_sun    = 8.0      # Bovy et al. (2009)
#R_sun    = 8.0      # Robin et al. 
R_sun   = 8.3      # Brunthaler et al. (2011)
#____________________________________________________________________________________________________
#Vclsr   =  220.    # Bovy et al. (2009)
Vclsr    = 239.0    # Brunthaler et al. (2011)   
Ulsr, Vlsr, Wlsr =  0.,   Vclsr,   0.
#____________________________________________________________________________________________________
#Usun, Vsun, Wsun = 10.4, 14.8  , 7.3                # Mihalas and Routhly (1968)
#Usun, Vsun, Wsun = 10.3, 20.5  , 7.5                # Erickson (1975) 
#Usun, Vsun, Wsun =  9. , 12.   , 7.                 # Mihalas and Binney  (1981)
#Usun, Vsun, Wsun = 10.3, 15.0  , 7.5                # Boesgaard and Tripico (1986)
#Usun, Vsun, Wsun = 11.0, 14.0  , 7.5                # Ratnatunga et al. (1989)
#Usun, Vsun, Wsun = 10.0, 5.23  , 7.17               # Dehnen & Binney (1998) and Hogg et al. (2005)
Usun, Vsun, Wsun = 11.1, 12.24 , 7.25               # Brunthaler et al. (2011) and Schonrich et al. (2010)
#Usun, Vsun, Wsun = 11.0, 12.0  , 7.0                # Robin et al. 
#____________________________________________________________________________________________________
Theta            = 122.932*factor                    # EPOCH 2000.
RAG, DECG        = 192.85948*factor, 27.12825*factor # EPOCH 2000.
#Theta           = 123.0*factor                      # EPOCH 1958.
#RAG, DECG       = 192.25*factor,  27.4*factor       # EPOCH 1958.
#____________________________________________________________________________________________________

output = open('RAVE_Kinematics_D83.out','a')
#output = open('Simu_m1407_selHRV-9-12.dat1_D83.out','a')
output.write('# RA_  DEC_  X  Y  Z  R_gal U V W VR  Vphi  Vb_star  GRV \n')

for i in np.arange(len(file_[:,0])):

# BGM

#	Vradial     = float(file_[i,7])
#	D           = float(file_[i,26])
#	l_,b_       = float(file_[i,22]), float(file_[i,23])
#	RA_,DEC_    = float(file_[i,24]), float(file_[i,25])
#	pmRA, pmDEC = float(file_[i,5]), float(file_[i,6])

# REAL DATA 
#
	Vradial     =  float(file_[i,6])
	D           =  float(file_[i,109])
	l_,b_       =  float(file_[i,4]), float(file_[i,5])
	RA_,DEC_    =  float(file_[i,2]), float(file_[i,3])
	pmRA, pmDEC =  float(file_[i,41]), float(file_[i,43])
	
	l, b        = l_*factor, b_*factor
	RA, DEC     = RA_*factor, DEC_*factor
	C1          = np.sin(DECG)*np.cos(DEC)-np.cos(DECG)*np.sin(DEC)*np.cos(RA-RAG)
	C2          = np.cos(DECG)*np.sin(RA-RAG)
	C           = np.matrix([[C1,C2],[-C2,C1]])
	pm          = np.matrix([[pmRA],[pmDEC]])
	pmG         = (1./np.cos(b))*C*pm
	
# Calculating Proper motions (mu_l, mu_b) and Velocities (Vl, Vb) in galactic coordinates from Radoslaw Polesky (2013)

	pml, pmb = np.cos(b)*pmG[0,0], pmG[1,0]
	Vl, Vb = D*k*np.cos(b)*pmG[0,0], D*k*pmG[1,0]

# Calculating X, Y and Z from Bond et al. (2010), ApJ, 716:1-29

	X     = R_sun-D*np.cos(l)*np.cos(b)
	Y     = -D*np.sin(l)*np.cos(b)
	Z     = D*np.sin(b)
	R_gal = np.sqrt((X*X)+(Y*Y)+(Z*Z)) # 3-dimension
	R     = np.sqrt((X*X)+(Y*Y))       # 2-dimension 

# Calculating Uobs, Vobs and Wobs from Bond et al. (2010), ApJ, 716:1-29

	Uobs = -Vradial*np.cos(l)*np.cos(b) + Vb*np.cos(l)*np.sin(b) + Vl*np.sin(l)
	Vobs = -Vradial*np.sin(l)*np.cos(b) + Vb*np.sin(l)*np.sin(b) - Vl*np.cos(l)
	Wobs =  Vradial*np.sin(b)           + Vb*np.cos(b)

# Calculating VR, Vphi ; Bond et al. (2010) 
# X axis is oriented toward l=180 degree, and the Y axis is oriented toward l = 270 degree, the disk rotates toward l = 90 degree.
	
	U   = Uobs - Usun - Ulsr
	V   = Vobs - Vsun - Vlsr
	W   = Wobs + Wsun - Wlsr

#	print U, V, W

#	U, V, W = -U, -V, W # Binney system

# Note that, in our adopted system, the disk has a prograde rotation Vphi = -220., retrograde rotation is indicated by Vphi > 0. Stars with VR > 0 move away from the Galactic center, and stars with Vz > 0 move toward the North Galactic Pole.	

	VX, VY, VZ = U, V, W
	VR         = (VX*(X/R))+(VY*(Y/R))
	Vphi       = (-VX*(Y/R))+(VY*(X/R))

# ______________________________________________________________________________________________________________
# Johnson and Soderblom (1987)
# We will use a right-handed coordinate system for U, V and W, so that they are positive in the directions of the Galactic center, Galactic rotation, and the North Galactic Pole (NGP), respectively.

	MatrixT    = sc.dot(sc.array([[sc.cos(Theta),sc.sin(Theta),0.],[sc.sin(Theta),-sc.cos(Theta),0.],[0.,0.,1.]]),sc.dot(sc.array([[-sc.sin(DECG),0.,sc.cos(DECG)],[0.,-1.,0.],[sc.cos(DECG),0.,sc.sin(DECG)]]),sc.array([[sc.cos(RAG),sc.sin(RAG),0.],[+sc.sin(RAG),-sc.cos(RAG),0.],[0.,0.,1.]])))
	MatrixA    = sc.array([[sc.cos(RA)*sc.cos(DEC),-sc.sin(RA),-sc.cos(RA)*sc.sin(DEC)],[sc.sin(RA)*sc.cos(DEC),sc.cos(RA),-sc.sin(RA)*sc.sin(DEC)],[sc.sin(DEC),0.,sc.cos(DEC)]])
	MatrixB    = sc.dot(MatrixT, MatrixA)
	VrVraVdec  = sc.array([Vradial, k*D*pmRA*np.cos(DEC), k*D*pmDEC])
	UVW        = sc.dot(MatrixB,VrVraVdec) 
	UU, VV, WW = UVW[0], UVW[1], UVW[2]
	UUc        = UU + Usun - Ulsr
	VVc        = VV + Vsun - Vlsr
	WWc        = WW + Wsun - Wlsr

# Calculating VLSR

	GRV        = Vradial + np.abs(Usun)*np.cos(l)*np.cos(b) + (np.abs(Vsun) + Vclsr)*np.sin(l)*np.cos(b) + np.abs(Wsun)*np.sin(b)
	Vb         = GRV / np.cos(b)
	
	output.write(str(RA_)+' '+str(DEC_)+' '+str(X)+' '+str(Y)+' '+str(Z)+' '+str(R_gal)+' '+str(UUc)+' '+str(VVc)+' '+str(WWc)+' '+str(VR)+' '+str(Vphi)+' '+str(Vb)+' '+str(GRV)+' \n')

output.close()
