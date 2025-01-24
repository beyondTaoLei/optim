#!/usr/bin/env python3
import os
import numpy as np
import matplotlib.pylab as plt

#Input paras
fdir='STD'
foptim=os.path.join(fdir,'optim.json')
optim   =eval(open(foptim).read())
niter_max =optim['niter_max']

#plot
iterno=np.arange(1,niter_max+1)
steps=np.zeros([niter_max,3],np.float32)
misfits=np.zeros(niter_max,np.float32)
for iter1 in range(niter_max):
    step   =eval(open(os.path.join(fdir,'iter'+str(iter1+1), 'step.json')).read())
    steps[iter1,:]=step['alpha1'],step['alpha2'],step['optim']
    misfits[iter1]=step['misfit0']

np.savetxt(os.path.join(fdir,'step_lens.dat'),steps,fmt='%.4e')
np.savetxt(os.path.join(fdir,'misfits.dat'),steps,fmt='%.4e')

#step length estimation
fig, ax = plt.subplots()
plt.plot(iterno,steps[:,2],'g',label='optim step',linewidth=2)
plt.plot(iterno,steps[:,0],'r.',label='1st step')
plt.plot(iterno, steps[:,1],'b.',label='2nd step')
plt.grid(True)
ax.set_xlim(0, niter_max)
ax.set_title('dynamatic step length estimation', fontsize=17)
plt.xlabel('Iterations', fontsize=15)
plt.ylabel('relative step length', fontsize=15)
handles, labels = plt.gca().get_legend_handles_labels()
plt.legend(handles,labels)
fig.savefig('fig/'+fdir+'_sl.eps',format='eps')
plt.show()

#misfit
fig, ax = plt.subplots()
plt.plot(iterno,misfits,'g-+',linewidth=2)
plt.grid(True)
ax.set_xlim(0, niter_max)
ax.set_title('misfit evolution', fontsize=17)
plt.xlabel('Iterations', fontsize=15)
plt.ylabel('misfit', fontsize=15)
fig.savefig('fig/'+fdir+'_misfit.eps',format='eps')
plt.show()

print("Finished...", __file__)
