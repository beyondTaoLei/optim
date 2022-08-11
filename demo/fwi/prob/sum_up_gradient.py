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
n3          =optim['n3']
n2          =optim['n2']
n1          =optim['n1']
fgrad       =optim['fgrad']
fout        =fgrad+'_0'

# interpretation
evn_info=pd.read_csv(fevn_info)
events=list(evn_info['station'])

# sum up
num = n1*n2*n3
sum1=np.zeros(num,np.float32)
for evn in events:
    fname=os.path.join(fdio,evn,'gradvp_shot1.bin')
    print(fname)
    if os.path.exists(fname):
        grad=np.fromfile(fname, np.float32)
        sum1+=grad
    else:
        print("no file: "+fname)

sum1.tofile(fout)
max1=np.abs(sum1).max()/10.

print("\nximage  n1=%d cmap=hsv2 bclip=%f wclip=%f <%s\n"% \
     (n1, max1, -max1, fout))
print("*"*60+"\n*Output file: "+fout)
print("Finished...", __file__)
