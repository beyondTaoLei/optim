#!/usr/bin/env python3
"""
Usage:
    python *.py
"""
import os
import sys
import numpy as np
import pandas as pd

#Input paras
fdir        =sys.argv[1]
foptim      =os.path.join(fdir,'optim.json')
optim       =eval(open(foptim).read())
fdio        =optim['fdio']
fevn_info   =optim['fevn_info']
nx          =optim['NX']
ny          =optim['NY']
idx         =optim['IDX']
idy         =optim['IDY']
fgrad       =optim['fgrad']
fout        =fgrad+'_0'

evn_info=pd.read_csv(fevn_info)
events=list(evn_info['station'])
n2=len(np.zeros(nx)[::idx])
n1=len(np.zeros(ny)[::idy])
sum1=np.zeros([n2,n1],np.float32)
for evn in events:
    fname=os.path.join(fdio,evn,'gradvp_shot1.bin')
    print(fname)
    if os.path.exists(fname):
        grad=np.fromfile(fname, np.float32).reshape([n2,n1])
        sum1+=grad
    else:
        print("no file: "+fname)

sum1.tofile(fout)
max1=np.abs(sum1).max()/10.

print("\nximage  n1=%d n2=%d cmap=hsv2 bclip=%f wclip=%f <%s\n"% \
     (n1, n2, max1, -max1, fout))
print("*"*60+"\n*Output file: "+fout)
print("Finished...", __file__)
