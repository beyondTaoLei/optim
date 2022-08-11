#!/usr/bin/env python3
"""
Usage:
    python *.py
"""
import os
import sys
import numpy as np
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter

#Input paras
fdir        =sys.argv[1]
foptim      =os.path.join(fdir,'optim.json')
optim       =eval(open(foptim).read())
n3g         =optim['n3g']
n2g         =optim['n2g']
n1g         =optim['n1g']
inc3        =optim['inc3']
inc2        =optim['inc2']
inc1        =optim['inc1']
n3          =optim['n3']
n2          =optim['n2']
n1          =optim['n1']
#fw          =optim['fw']
#free_surf   =optim['free_surf']
smooth_size =optim['smooth_size'] #5
taper       =optim['taper']
ftaper      =optim['ftaper']
fgrad       =optim['fgrad']
fin         =fgrad+'_0'
fout        =fgrad+'_1'

#a. interpolate
grad=np.fromfile(fin, np.float32)
grid3f, grid2f, grid1f = np.mgrid[0:(n3g-1):n3g*1j, 0:(n2g-1):n2g*1j, 0:(n1g-1):n1g*1j]
grid3i, grid2i, grid1i = np.mgrid[0:(n3-1):n3*1j, 0:(n2-1):n2*1j, 0:(n1-1):n1*1j]
grid3i, grid2i, grid1i = inc3*grid3i, inc2*grid2i, inc1*grid1i
points  = (grid3i.flatten(), grid2i.flatten(), grid1i.flatten())
pointsN = (grid3f.flatten(), grid2f.flatten(), grid1f.flatten())
grad = griddata(points, grad, pointsN, method='nearest')

#efficitive grids due to boundary saving strategy
#here we extend 1 grid on the left and right sides
grad = grad.reshape([n3g, n2g, n1g])
#if fw > 0:
    ## n3
    #if n3g > 1:
        #grad[:fw+1,:,:]=0.0
        #grad[n3g-fw-1:,:,:]=0.0
    ## n2
    #grad[:,:fw+1,:]=0.0
    #grad[:,n2g-fw-1:,:]=0.0
    ## n1
    #grad[:,:,n1g-fw-1:]=0.0
    #if free_surf==1:
        #grad[:,:,:2]=0.0
    #else:
        #grad[:,:,:fw+1]=0.0

#b. smooth
print("smooth...")
sigmas=smooth_size/np.sqrt(8)
grad = gaussian_filter(grad, sigmas, mode='nearest')#sigma

#c. taper
if taper==1:
    print("taper...")
    if not os.path.exists(ftaper):
        raise ValueError("no taper file "+ftaper)
    taper=np.fromfile(ftaper, np.float32)
    grad=grad.flatten()*taper

grad.astype(np.float32).tofile(fout)
max1=np.abs(grad).max()/10.

print("\nbefore processing...")
print("\nximage  n1=%d n2=%d cmap=hsv2 bclip=%f wclip=%f <%s\n"% \
     (n1, n2, max1, -max1, fin))
print("after processing...")
print("\nximage  n1=%d n2=%d cmap=hsv2 bclip=%f wclip=%f <%s\n"% \
     (n1g, n2g, max1, -max1, fout))
print("*"*60+"\n*Output file: "+fout)
print("Finished...", __file__)
