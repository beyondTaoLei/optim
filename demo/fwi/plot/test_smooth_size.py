#!/usr/bin/env python3
import os
import sys
import numpy as np
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter
import matplotlib.pylab as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

#Input paras
nxr         =int(sys.argv[1])
nzr         =int(sys.argv[2])
grad_filt   =int(sys.argv[3]) #[-3*sigma, 3*sigma] will be covered, sigma = (1_of_4_wavelength_in_grids-1)/2.0/3.0
oxr, ozr = 0.0, 0.0
dxr, dzr = 1.0, 1.0
# nxr, nzr = 1000, 300
# grad_filt = 20 #sigma
px = nxr//2
pz = nzr//2

data = np.zeros([nxr, nzr],np.float32) # bin grid starting from left bottom corner
# data[px,pz-50]= 1.0
# data[px+50,pz]= 1.0
data[px,pz]= 1.0
datasm = gaussian_filter(data, sigma=grad_filt)

#update models
max_desc =np.amax(np.fabs(datasm))
vmin = 0.0
vmax = 0.2*max_desc

# line maker
datasm[:,pz-grad_filt]= 100.0
datasm[:,pz+grad_filt]= 100.0
datasm[px-grad_filt,:]= 100.0
datasm[px+grad_filt,:]= 100.0
datasm[:,pz-3*grad_filt]= 100.0
datasm[:,pz+3*grad_filt]= 100.0
datasm[px-3*grad_filt,:]= 100.0
datasm[px+3*grad_filt,:]= 100.0
datasm[:,pz]= 0.1*max_desc
datasm[px,:]= 0.1*max_desc

# show figure or not
flag_show    = 1
# figure name to save
fignm   = 'sm.png'
# figure size to save
figsize = [6,3]
# figure resolution to save
figdpi  = 150
name = 'grad_filt'

##########################################
X=np.arange(nxr)*dxr+oxr
Z=np.arange(nzr)*dzr+ozr
coord0, coord1 = np.meshgrid(X, Z, indexing='ij')
# coord0 = coord0/1e3
# coord1 = coord1/1e3
print('axis 0:', coord0.min(),coord0.max())
print('axis 1:', coord1.min(),coord1.max())

# media show
fig = plt.figure(dpi=figdpi,figsize=(figsize[0],figsize[1]))#,figsize=(figsize[0],figsize[1]))
ax = plt.gca()
im = plt.pcolor(coord0, coord1, datasm, shading='auto', cmap = 'jet', vmin=vmin, vmax=vmax)
plt.xlabel('X (grid)')
plt.ylabel('Z (grid)')
plt.axis('image')
plt.gca().invert_yaxis()
plt.title(name)
#legend
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
clb = plt.colorbar(im, cax=cax)
plt.savefig(fignm)
if flag_show:
  plt.show()

print("Finished...", __file__)
