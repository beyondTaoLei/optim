#!/usr/bin/env python3
import os
import sys
import numpy as np  
import pandas as pd
import matplotlib.pylab as plt

#Input paras
fdir        =sys.argv[1]
fhist = os.path.join(fdir, 'history.csv')
df = pd.read_csv(fhist)
iters = df['iter']
misfits = df['misfit']

#figure
params = {"figure.figsize":(6,4),
          'axes.labelsize': 9,
          'axes.titlesize': 9,
          'xtick.labelsize': 9,
          'ytick.labelsize': 9,
          'legend.fontsize': 7
          }
plt.rcParams.update(params)

# media show
fig, ax1 = plt.subplots(dpi=150)
im1 = ax1.plot(iters, misfits, 'k', marker='o')
ax1.set_title('misfit evolution')
ax1.set_xlabel('Iteration')
ax1.set_ylabel('Misfit')
# plt.xticks(np.arange(1, 11, 1))
# plt.yticks(np.arange(0.5, 6, 0.5)*1e6)
ax1.grid(True)


plt.savefig('fig/misfit.png')
plt.show()
    
print("Finished...", __file__)    
    
        
