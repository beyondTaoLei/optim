#!/usr/bin/env python3
"""
Usage:
    python *.py
"""
import os
import numpy as np
from util import medfilt2d_enhanced
from scipy.interpolate import RectBivariateSpline

#Input paras
fdir        =sys.argv[1]
foptim      =os.path.join(fdir,'optim.json')
optim       =eval(open(foptim).read())
n2g         =optim['NX']
n1g         =optim['NY']
FW          =optim['FW']
FREE_SURF   =optim['FREE_SURF']
id2         =optim['IDX']
id1         =optim['IDY']
smooth_size =optim['smooth_size'] #5
taper       =optim['taper']
ftaper      =optim['ftaper']
fgrad       =optim['fgrad']
fin         =fgrad+'_0'
fout        =fgrad+'_1'

#a. interpolate
n2=len(np.zeros(n2g)[::id2])
n1=len(np.zeros(n1g)[::id1])
grad0=np.fromfile(fin, np.float32).reshape([n2,n1])
if id2!=1 or id1!=1:
    print("resample...")
    #interpolate2d from inversion grid to fd grid
    grd2_den, grd1_den=np.arange(n2g), np.arange(n1g)
    grd2_sam, grd1_sam=np.arange(n2g)[::id2], np.arange(n1g)[::id1]
    fun = RectBivariateSpline(grd2_sam, grd1_sam, grad0, kx=1, ky=1)
    grad0 =fun(grd2_den, grd1_den)
else:
    print("resample... no")
    grad0=grad0

#efficitive grids due to boundary saving strategy
#here we extend 1 grid on the left and right sides
grad0[:FW+1,:]=0.0
grad0[n2g-FW-1:,:]=0.0
if FREE_SURF==1:
    grad0[:,:2]=0.0
else:
    grad0[:,:FW+1]=0.0
grad=grad0
#b. smooth
print("smooth...")
grad=medfilt2d_enhanced(grad, smooth_size)
#c. taper
if TAPER==1:
    print("taper...")
    if not os.path.exists(ftaper):
        raise ValueError("no taper file "+ftaper)
    taper=np.fromfile(ftaper, np.float32).reshape([n2g,n1g])
    grad=grad*taper
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
