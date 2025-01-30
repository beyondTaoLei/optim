#!/usr/bin/env python3
"""
Usage:
    python *.py
"""
import os
import sys
import json
import numpy as np
from scipy.ndimage import gaussian_filter

#Input paras
fdir        =sys.argv[1]
foptim      =os.path.join(fdir,'optim.json')
optim       =eval(open(foptim).read())
nx, nz      =optim['nxzfd']
grad_filt   =optim['grad_filt'] #5
grad_taper  =optim['grad_taper']
ftaper      =optim['grad_taper_file']
fgrad       =optim['fgrad']
eps_scale   =optim['eps_scale']
fmodfd      =optim['fmodfd']
fin         =fgrad+'_0'

#a. interpolate
grad = np.fromfile(fin, np.float32).reshape([nx, nz])

#b. smooth
print("smooth...")
grad = gaussian_filter(grad, grad_filt, mode='nearest')#sigma

#c. taper
if grad_taper == 1:
  print("taper...")
  if not os.path.exists(ftaper):
      raise ValueError("no taper file "+ftaper)
  taper = np.fromfile(ftaper, np.float32).reshape([nx, nz])
  grad *= taper

#d. scaling
if 'SCALING_FACTOR' not in optim:
    mod =np.fromfile(fmodfd, np.float32)
    mod_max=np.abs(mod).max().item()
    grad_max=np.abs(grad).max().item()
    SCALING_FACTOR=eps_scale*mod_max/grad_max
    SCALING_FACTOR=1.0/SCALING_FACTOR
    print("*** generate SCALING_FACTOR at first iteration ***")
    print("SCALING_FACTOR= %e"%SCALING_FACTOR)
    print("eps_scale=%f,mod_max=%.2f,grad_max=%f\n"%(eps_scale,mod_max,grad_max))
    optim.update({'modmax0':   mod_max, \
                  'gradmax0':  grad_max,\
                  'SCALING_FACTOR':SCALING_FACTOR
            })
    fout2 = open(foptim,'w')
    fout2.write(json.dumps(optim,indent=4))
    fout2.close()
else: #read
    SCALING_FACTOR =optim['SCALING_FACTOR']
    print("*** read SCALING_FACTOR from file: "+foptim+" ***")
    print("SCALING_FACTOR= %e"%SCALING_FACTOR)

max0=np.abs(grad).max()/10.
grad=grad/SCALING_FACTOR
max1=np.abs(grad).max()/10.
grad.astype(np.float32).tofile(fgrad)

print("\nbefore processing...")
print("\nximage  n1=%d cmap=hsv2 bclip=%f wclip=%f <%s\n"% \
     (nz, max0, -max0, fin))
print("after processing...")
print("\nximage  n1=%d cmap=hsv2 bclip=%f wclip=%f <%s\n"% \
     (nz, max1, -max1, fgrad))
print("*"*60+"\n*Output file: "+fgrad)
print("Finished...", __file__)