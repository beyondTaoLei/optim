#!/usr/bin/env python3
import os
import sys
import shutil
import numpy as np
import pandas as pd
from mpi4py import MPI
from disba import PhaseDispersion
from scipy.ndimage.filters import gaussian_filter
from util import BLOCK_LOW, BLOCK_HIGH, BLOCK_SIZE

#initialize MPI
comm=MPI.COMM_WORLD
NP = comm.Get_size()
MYID = comm.Get_rank()

#input paras
fdir        =sys.argv[1]
foptim      =os.path.join(fdir,'optim.json')
optim       =eval(open(foptim).read())
fmod        =optim['fmodfd']
fperiods    =optim['fperiods']
fcmp        =optim['fcmp']
wdir        =optim['wdir']
algo        =optim['algo'] #'dunkin'
wavetype    =optim['wavetype']#'love' or 'rayleigh'
nper        =optim['nper']
ntraces     =optim['ntraces']
n1          =optim['n1g']
d1          =optim['d1g']
maxmode     =optim['maxmode']
modes       =maxmode + 1

# vs perturbation (km/s): avoid the homogeneous media
if 'pertb' in optim.keys():
    pertb = optim['pertb']
else:
    pertb = 0.001
pertbs = (-1.0)**np.arange(n1)*pertb

# make directory
if MYID == 0:
    if os.path.exists(wdir):
        shutil.rmtree(wdir, ignore_errors=True)
    os.makedirs(wdir,exist_ok=True)
    os.makedirs(os.path.join(wdir,'part'),exist_ok=True)
comm.Barrier()

# interpretation
df = pd.read_csv(fperiods)
periods = np.array(df['T0'])

# read indexes of CMP
df = pd.read_csv(fcmp)
offsets = df['offset']
idx1  =BLOCK_LOW(MYID,NP,ntraces)
idx2  =BLOCK_HIGH(MYID,NP,ntraces)
mynum =BLOCK_SIZE(MYID,NP,ntraces)

# split models
modloc=np.zeros([3, mynum, n1], np.float32) #period from short to long
fvp = fmod[:-3]+'.vp'
fvs = fmod[:-3]+'.vs'
frho = fmod[:-3]+'.rho'
for idx in range(idx1, idx2+1):
    i = idx - idx1
    offset1 = offsets[idx]*n1*4
    vp1d = np.fromfile(fvp, np.float32, n1, offset=offset1)
    vs1d = np.fromfile(fvs, np.float32, n1, offset=offset1)
    rho1d = np.fromfile(frho, np.float32, n1, offset=offset1)
    modloc[0,i,:] = vp1d/1000.0
    modloc[1,i,:] = vs1d/1000.0 + pertbs
    modloc[2,i,:] = rho1d/1000.0

# save local models
modloc.tofile(os.path.join(wdir, 'part', 'mod_id%02d.bin'%MYID))

# run FD
model=np.zeros([n1, 4])
coordz1d = np.arange(n1)*d1
z = coordz1d/1000.0 #'kilometers'
tn = np.append(np.diff(z), [0, ])
model[:,0] = tn
phvel=np.zeros([modes, mynum, nper], np.float32) #period from short to long
for idx in range(idx1, idx2+1):
    i = idx - idx1
    model[:,1] = modloc[0,i,:]
    model[:,2] = modloc[1,i,:]
    model[:,3] = modloc[2,i,:]
    pd0 = PhaseDispersion(*model.T, algorithm=algo)
    for mode in range(modes):
        try:
            cp = pd0(periods, mode, wave=wavetype)
            pers, cs = cp.period, cp.velocity
            nper0 = len(pers)
        except:
            #smooth model
            print("smoothing at index: ", offsets[idx])
            model[:,2] = gaussian_filter(model[:,2], 6, mode='nearest')
            pd0 = PhaseDispersion(*model.T, algorithm=algo)
            #retry
            try:
                cp = pd0(periods, mode, wave=wavetype)
                pers, cs = cp.period, cp.velocity
                nper0 = len(pers)
            except:
                nper0 = 0
                print("Warning: an exception occurred at index: ", offsets[idx])
        #if nper0 > 0:
        phvel[mode, i, :nper0] = cs

phvel[:] *= 1000.0 #m/s
for mode in range(modes):
    fnm = os.path.join(wdir, 'part', 'vmap_M%1d_id%02d.bin'%(mode, MYID))
    phvel[mode,:,:].tofile(fnm)

comm.Barrier()
#gather phvel
for mode in range(modes):
    if mode%NP==MYID:
        fout = os.path.join(wdir, 'vmap_M%1d.bin'%mode)
        with open(fout, 'wb') as fp:
            for id1 in range(NP):
                fnm = os.path.join(wdir, 'part', 'vmap_M%1d_id%02d.bin'%(mode, id1))
                part=np.fromfile(fnm, np.float32)
                part.tofile(fp)

if MYID == 0:
    print("Finished...", __file__)