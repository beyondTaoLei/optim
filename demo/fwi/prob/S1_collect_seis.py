#!/usr/bin/env python3
"""
collect seismograms from all sources
"""
import os
import sys
import h5py
import obspy
import numpy as np
import pandas as pd
from mpi4py import MPI
from shutil import copyfile

#initialize MPI
comm=MPI.COMM_WORLD
NP = comm.Get_size()
MYID = comm.Get_rank()

#Input paras
fdir        =sys.argv[1]
flag        =sys.argv[2] #[all/test]: all shots or test shots
foptim      =os.path.join(fdir,'optim.json')
optim       =eval(open(foptim).read())
delta       =optim['delta']
nt          =optim['nt']
comp        =optim['comp'] #'Vz'
fevn_info   =optim['fevn_info']
fdio        =optim['fdio']
fdata       =optim['fdata']

#read stn info.
df = pd.read_csv(fevn_info)
if flag == 'test':
  cond = df['test'] == 1
  df = df.loc[cond]

test = df['test'].to_numpy()
evns = list(df['station'])
Scoord = np.vstack([df['x'], df['z']]).transpose() # meter
nsrc = len(evns)
comm.Barrier()

#output file
fout = os.path.join(fdata, 'syn0.h5')
hf = h5py.File(fout, 'w', driver='mpio', comm=MPI.COMM_WORLD)
# collect seismograms
for isrc in range(nsrc):
  kevnm = evns[isrc]
  fnm = os.path.join(fdio, kevnm, 'receiver.dat')
  Rcoord = np.loadtxt(fnm)
  ntr, tmp = Rcoord.shape
  Rxyz = np.vstack([Rcoord[:,0], Rcoord[:,0]*0.0, Rcoord[:,1]]).transpose()
  Sxyz = np.array([Scoord[isrc, 0], 0.0, Scoord[isrc, 1]])
  #write
  grp = hf.create_group(kevnm)
  grp.attrs['test'] = test[isrc]
  grp.attrs['b'] = 0.0
  grp.attrs['delta'] = delta
  grp.attrs['nt'] = nt 
  grp.attrs['ntr'] = ntr
  grp.attrs['comp'] = comp
  dset = grp.create_dataset('S', data=Sxyz, dtype='f') # 3
  dset = grp.create_dataset('R', data=Rxyz, dtype='f') # [ntr, 3]
  ## require space
  for cmp in comp:
    dset = grp.require_dataset(cmp, shape=(ntr, nt), dtype='f')

comm.Barrier()

# collect seismograms
for isrc in range(MYID, nsrc, NP):
  kevnm = evns[isrc]
  grp = hf[kevnm]
  ntr   = grp.attrs['ntr']
  nt    = grp.attrs['nt']
  for cmp in comp:
    data = np.zeros([ntr, nt], np.float32)
    fnm = os.path.join(fdio, kevnm, 'seis_'+cmp.lower()+'_shot1.su') #['p','vx','vy']
    if not os.path.exists(fnm):
      print(fnm + ": No data")
    else:
      print(fnm)
      seis = obspy.read(fnm)
      for itr in range(ntr):
        data[itr,:] = seis[itr].data
    #write
    dset = grp[cmp]
    dset[:] = data

hf.close()

#back up the synthetic data at current model
if MYID == 0 and flag != 'test':
  copyfile(fout, fout.replace('syn0.h5','syn_bk.h5'))

print("Finished...", __file__)
