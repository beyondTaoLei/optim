#!/usr/bin/env python3
import os
import sys
import numpy as np
import pandas as pd
from mpi4py import MPI
sys.path.append('/home/tao/seis_software/ATT_surfwaves/forward')
from dpcxx import GradientEval

#initialize MPI
comm=MPI.COMM_WORLD
NP = comm.Get_size()
MYID = comm.Get_rank()

#input paras
fdir        =sys.argv[1]
mode        =int(sys.argv[2]) #0, 1
foptim      =os.path.join(fdir,'optim.json')
optim       =eval(open(foptim).read())
fmod        =optim['fmodfd']
fperiods    =optim['fperiods']
wdir        =optim['wdir']
wavetype    =optim['wavetype']
nper        =optim['nper']
ntraces     =optim['ntraces']
n1          =optim['n1g']
d1          =optim['d1g']
lastele     =3

# interpretation
df = pd.read_csv(fperiods)
periods = np.array(df['T0'])
freqs = 1.0/periods

# read phase velocity
fnm = os.path.join(wdir, 'part', 'vmap_M%1d_id%02d.bin'%(mode, MYID))
fsize=os.path.getsize(fnm)
if fsize%(4*nper)!=0: #sizeof(float)
    raise ValueError("Check the size of file: " + fnm)
mynum = fsize//(4*nper)
phvel = np.fromfile(fnm, np.float32).reshape([mynum, nper])
phvel[:] /= 1000.0 #km/s

# read model
fnm = os.path.join(wdir, 'part', 'mod_id%02d.bin'%MYID)
modloc = np.fromfile(fnm, np.float32).reshape([3, mynum, n1])
model = np.zeros([n1, 5])
model[:,0] = np.arange(n1) + 1.0
coordz1d = np.arange(n1)*d1
model[:,1] = coordz1d/1000.0 #'kilometers'

# run kernel generation
kernels = np.zeros([nper, mynum, n1], np.float32)
for i in range(mynum):
    model[:,2] = modloc[2,i,:] #rho
    model[:,3] = modloc[1,i,:] #vs
    model[:,4] = modloc[0,i,:] #vp
    cs = phvel[i, :]
    gradEval = GradientEval(model, wavetype)
    for f1, c1 in zip(freqs, cs):
        if c1 > 0.0:
            iper=(np.abs(periods - 1.0/f1)).argmin()
            kernels[iper, i, :] = gradEval.compute(f1, c1)
            # reset the last two elements
            kernels[iper, i, -lastele:] =kernels[iper, i, -lastele]

# save kernels
for iper in range(nper):
    fname='ker_M%1d_T%.2f_id%02d.bin'%(mode, periods[iper], MYID)
    fout=os.path.join(wdir, 'part', fname)
    kernels[iper, :, :].tofile(fout)

comm.Barrier()
#gather kernel
for iper in range(nper):
    if iper%NP==MYID:
        fname='ker_M%1d_T%.2f.bin'%(mode, periods[iper])
        fout=os.path.join(wdir, fname)
        with open(fout, 'wb') as fp:
            count=0
            for id1 in range(NP):
                fker0='ker_M%1d_T%.2f_id%02d.bin'%(mode, periods[iper], id1)
                fker=os.path.join(wdir, 'part', fker0)
                part=np.fromfile(fker, np.float32)
                where_are_NaNs=np.isnan(part)
                part[where_are_NaNs] = 0.0 #unknown error
                part.tofile(fp)
                count+=np.sum(where_are_NaNs)
            if count > 0:
                print('NaN percentage: %s %.6f'%(fname, count/(ntraces*n1)))

if MYID == 0:
    print("Finished...", __file__)