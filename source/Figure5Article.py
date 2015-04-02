#!/usr/bin/python


import numpy as np
import scipy as sc
import pylab as plt
import matplotlib as mpl

f=plt.figure(1)#, (30, 40))
#f.subplots_adjust(hspace=0.74,wspace=0.24, bottom = 0.11, top = 0.86, left = 0.12, right = 0.94)

#etiq = ['2','25','1','5']
#etiq = ['1','2','3','4']
Nlabel =20

for i in np.arange(89): 
	data = sc.genfromtxt(str(int(i)+1)+'.vt')
	d    = data[:,0]*1000. 
	v    = data[:,1]*10./1000.


	binx1 = int(1000./50.)
	biny1 = int(10.)

	ax = f.add_subplot(15,6,int(i)+1)
	hist,xedges,yedges = np.histogram2d(v,d,bins=(binx1,biny1))
	aspectratio = 1.0*(1000.)/(100.)
	masked = np.ma.masked_where(hist==0, hist)
#	cmap_multicolor = mpl.cm.gray_r.set_bad('w', 1)
	r_=ax.imshow(masked.T,extent=[0.,1.,0.,117.],interpolation='nearest',aspect='auto',cmap=mpl.cm.gray_r,origin='lower')
	ae_=plt.colorbar(r_,shrink=0.85)
	#ax.text(700.,5.,'#'+str(etiq[i]),color='red',fontsize=30)
	ae_.ax.tick_params(labelsize = Nlabel)
	ax.set_xlabel(r'$V_{rel}$ (km/s)',fontsize=Nlabel)
	ax.set_ylabel(r'$d_{min}$ (pc)',fontsize=Nlabel)
	ax.tick_params(labelsize = Nlabel)
	ax.set_xlim(0.,1.)
	ax.set_ylim(0.,117.)
	ax.text(850./1000., 20., int(i)+1, color ="black",fontsize=17)
#	ax.contour(masked.T/1.E5*100., 3, colors='k', linewidths=1,extent=[np.min(v),np.max(v),np.min(d),np.max(d)],linestyles='dashed')
#	ax.axhline(y=86.23,ls='-' ,lw=2,color='red')  # Harris et al. (1999) and Da Costa et al. (2008), 5.2 kpc
#	ax.axhline(y=107.87,ls='-',lw=2,color='cyan') # Fernandez-Trincado et al. (2015), 5.15 kpc , 1.2deg Marconi et al. (2014)
#	ax.axhline(y=115.20,ls='-',lw=2,color='green') # Del Principe et al. (2006), 5.5 kpc, 1.2deg Marconi et al. (2014)
#	ax.axhline(y=116.67,ls='--',lw=2,color='red') # Navarrete et al. (2015), 5.57 kpc, 1.2deg Marconi et al. (2014)


plt.show()
