#!/usr/bin/env python3
"""
Usage:
    python *.py
"""
import os
import sys
import json
import numpy as np

#Input paras
fdir        =sys.argv[1]
foptim      =os.path.join(fdir,'optim.json')
optim       =eval(open(foptim).read())
eps_scale   =optim['eps_scale']
fmod        =optim['fmodfd']
fgrad       =optim['fgrad']
n1g         =optim['n1g']
fin         =fgrad+'_1'

grad=np.fromfile(fin, np.float32)
if 'SCALING_FACTOR' not in optim:
    mod =np.fromfile(fmod, np.float32)
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

grad=grad/SCALING_FACTOR
grad.tofile(fgrad)
max1=np.abs(grad).max()/10.

print("\nximage  n1=%d cmap=hsv2 bclip=%f wclip=%f <%s\n"% \
     (n1g, max1, -max1, fgrad))
print("*"*60+"\n*Output file: "+fgrad)
print("Finished...", __file__)
