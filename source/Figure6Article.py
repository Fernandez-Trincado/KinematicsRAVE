#!/usr/bin/python


import numpy as np
import scipy as sc
import pylab as plt


Nstar  = 89.
f      = plt.figure(1)
Nlabel = 20

for i in np.arange(int(Nstar)):

	data = sc. genfromtxt('dat1.orbs.bprol.omcest'+str(int(i)+1)+'.marzo2015.ax.gauss.1')

	t    = (data[:,16]*(-9.78462E7))/(1.E9)
	rt   = data[:,17]*1000.
	vr   = data[:,18]*10.


	ax   = f.add_subplot(15,6,int(i)+1)
	
	maskrt = (rt <= 117.)	
	maskvr = (vr[maskrt] < 200.)


	hist_, bins_ = np.histogram(t[maskrt],bins=10,range=[0.,1.])
#	plt.hist(t[maskrt],lw=3., histtype='stepfilled',color='grey',edgecolor='grey',normed=True,bins=10, range=[0.,1.])	

	plt.plot(np.append([0.],np.linspace(0.,1.,10.)) , np.append([0.],hist_) )


#	plt.plot(bins_[:-1],hist_/float(np.max(hist_)))


#	bin_ = (np.max(hist_)/10.)
#	plt.xlabel('Time (Gyr)')
#	plt.ylabel('N')



#	plt.text(0.85, bin_, int(i)+1, color ="black",fontsize=17)	
#
#	if int(len(t[maskrt][maskvr])) > 0: 
#
#		plt.hist(t[maskrt][maskvr], lw=3., histtype='stepfilled', fill=False,edgecolor='black',normed=False,bins=10)
#
#	else: pass


plt.show()


