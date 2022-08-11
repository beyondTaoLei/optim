#!/usr/bin/env python3
"""
Usage:
    python *.py
"""
import os
import sys
import obspy
import pandas as pd
from mpi4py import MPI
from shutil import copyfile
import obspy.signal.filter as flt

#initialize MPI
comm=MPI.COMM_WORLD
NP = comm.Get_size()
MYID = comm.Get_rank()

#Input paras
fdir        =sys.argv[1]
flist       =sys.argv[2]
fc          =float(sys.argv[3])   #low cut-off frequency in Hz
foptim      =os.path.join(fdir,'optim.json')
optim       =eval(open(foptim).read())
order       =optim['order']#4
zero_phase  =optim['zero_phase']
if zero_phase==0:
    zero_phase = False  
else:
    zero_phase = True

# interpretation
df = pd.read_csv(flist)
fin = list(df['fin'])
fout = list(df['fout'])
if fc > 0:
    for i in range(MYID, len(fin), NP):
        if not os.path.exists(fin[i]):
            raise ValueError("no file: "+fin[i])
        
        print("Butterworth filtering(fc order= %.2f %d ) at %s"%(fc, order, fin[i]))
        seis=obspy.read(fin[i])
        delta=seis[0].stats.delta
        ntr=seis.count()
        #Apply a zero-phase lowpass to the trace
        for itr in range(ntr):
            seis[itr].data[:]=flt.lowpass(seis[itr].data,fc,1.0/delta, order, zerophase=zero_phase)
        seis.write(fout[i])
else:
    for i in range(MYID, len(fin), NP):
        print("copy file: from %30s to %30s"%(fin[i], fout[i]))
        copyfile(fin[i],  fout[i])
        

print("Finished...", __file__)
        
