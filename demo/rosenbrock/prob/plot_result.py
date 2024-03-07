#!/usr/bin/env python3
import os,sys
import numpy as np
from rosenbrock import func
import matplotlib.pylab as plt
import matplotlib.cm as cm

#Input paras
fdir        =sys.argv[1]#'optim_LBFGS'
niter_max   =int(sys.argv[2])#16#25#50#1000#optim['niter_max']

#plot
invs=np.zeros([niter_max,2],np.float32)
for iter1 in range(niter_max):
    invs[iter1, :]=np.fromfile(os.path.join(fdir,'iter'+str(iter1+1), 'mod.bin'), np.float32)
    
np.savetxt(os.path.join(fdir,'result.dat'),invs,fmt='%.4e')

#invs=np.loadtxt('optim/result.dat')
max1=1.1
X=np.arange(100)/100*max1
Y=np.arange(100)/100*max1
data=np.zeros([len(X),len(Y)])
for i in range(len(X)):
    for j in range(len(Y)):
        data[i,j]=func([X[i],Y[j]])

#plt.show()
cmap=cm.bwr
fig, ax = plt.subplots(figsize=(12,10))
data=data.transpose()#for ximage
im = ax.imshow(data, interpolation='bilinear', cmap=cmap, aspect='auto', vmax=8, vmin=0,extent=[0,max1,max1,0.])
ax.set_title('rosenbrock 2D-'+fdir[6:], fontsize=17)
plt.xlabel('X', fontsize=15)
plt.ylabel('Y', fontsize=15)
plt.plot(invs[:,0], invs[:,1],'w+-',linewidth=2)
ax.text(invs[-1,0], invs[-1,1],niter_max,c='k')
plt.plot([1.0], [1.0],'ro',label='true')
plt.plot([0.25], [0.25],'ko',label='initial')
handles, labels = plt.gca().get_legend_handles_labels()
plt.legend(handles,labels)
fig.colorbar(im)
fig.savefig('fig/'+fdir+'.eps',format='eps')
plt.show()
print("Finished...", __file__)
