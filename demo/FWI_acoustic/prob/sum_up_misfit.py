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
flag        =sys.argv[2] #[all/test]: all shots or test shots
foptim      =os.path.join(fdir,'optim.json')
optim       =eval(open(foptim).read())
fdio        =optim['fdio']
fevn_info   =optim['fevn_info']
lnorm       =optim['lnorm']
fcost       =optim['fcost']

evn_info=pd.read_csv(fevn_info)
test=np.array(evn_info['test'])
if flag=="all":
    print("sum up misfit from all shots")
    events=list(evn_info['station'])
else:
    print("sum up misfit from the test shots")
    print("and the number of test shots is", test.sum())
    events=list(evn_info['station'][test==1])

L2sum, energy_sum =0., 0.
for evn in events:
    fname=os.path.join(fdio,evn,"misfit.json")
    print(fname)
    if not os.path.exists(fname):
        raise ValueError("no file: "+fname)
    misf_paras = eval(open(fname).read())
    if lnorm==1:
        L2sum +=misf_paras['misfit']
        energy_sum +=misf_paras['energy']
    else:
        raise ValueError("we will implement other norm in future!")

if lnorm==1:
    result=L2sum/energy_sum
fp = open(fcost,'w')
fp.write('%.6e\n'% result)
fp.close()
print("*"*60+"\n*Output file: "+fcost)
print("Finished...", __file__)
