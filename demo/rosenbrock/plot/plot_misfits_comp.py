#!/usr/bin/env python3
import os
import numpy as np
import matplotlib.pylab as plt

#plot
misfits0=np.zeros(50,np.float32)
misfits1=np.zeros(50,np.float32)
misfits2=np.zeros(50,np.float32)
misfits3=np.zeros(50,np.float32)
for iter1 in range(50):
    step   =eval(open(os.path.join('optim_PSTD','iter'+str(iter1+1), 'step.json')).read())
    misfits0[iter1]=step['misfit0']

for iter1 in range(50):
    step   =eval(open(os.path.join('optim_PNLCG','iter'+str(iter1+1), 'step.json')).read())
    misfits1[iter1]=step['misfit0']
    step   =eval(open(os.path.join('optim_LBFGS','iter'+str(iter1+1), 'step.json')).read())
    misfits2[iter1]=step['misfit0']
    step   =eval(open(os.path.join('optim_TRN','iter'+str(iter1+1), 'step.json')).read())
    misfits3[iter1]=step['misfit0']

#misfit
fig, ax = plt.subplots()
plt.plot(np.arange(50)+1,misfits0[:50],'k',linewidth=2,label='steepest descent')
plt.plot(np.arange(50)+1,misfits1,'b',linewidth=2,label='nonlinear CG')
plt.plot(np.arange(50)+1,misfits2,'g',linewidth=2,label='l-BFGS')
plt.plot(np.arange(50)+1,misfits3,'r',linewidth=2,label='truncated Newton')
plt.grid(True)
ax.set_xlim(0, 50)
ax.set_title('misfit evolution', fontsize=17)
plt.xlabel('Iterations', fontsize=15)
plt.ylabel('misfit', fontsize=15)
handles, labels = plt.gca().get_legend_handles_labels()
plt.legend(handles,labels)
fig.savefig('fig/'+'misfit_all.eps',format='eps')
plt.show()

print("Finished...", __file__)
