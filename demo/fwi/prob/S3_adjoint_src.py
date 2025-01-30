#!/usr/bin/env python3
"""
calculate the adjoint source
"""
import os
import sys
import h5py
import numpy as np
import pandas as pd
from mpi4py import MPI
from shutil import copyfile
from util import write_su
from misfit_funs import misf_wfm, misf_tt

#initialize MPI
comm=MPI.COMM_WORLD
NP = comm.Get_size()
MYID = comm.Get_rank()

#Input paras
fdir          =sys.argv[1]
flag          =int(sys.argv[2])
foptim        =os.path.join(fdir,'optim.json')
optim         =eval(open(foptim).read())
fdio          =optim['fdio']
fdata         =optim['fdata']
Lnorm         =optim['Lnorm']
fmisfit       =optim['fmisfit']
fcost         =optim['fcost']
# fevn_info     =optim['fevn_info']

'''
flag:
1: calculate adjoint source and misfit function
2: same as case 1, but do not store the adjoint source, only for case of test=1
3: collect misfit from misfit.csv for step = 0
'''
fobs = os.path.join(fdata, 'obs.h5')
fsyn = os.path.join(fdata, 'syn.h5')
fadj = os.path.join(fdata, 'adj.h5')
hfobs = h5py.File(fobs, 'r', driver='mpio', comm=MPI.COMM_WORLD)
hfsyn = h5py.File(fsyn, 'r', driver='mpio', comm=MPI.COMM_WORLD)
evns = list(hfsyn.keys()) # should use syn info.
nsrc = len(evns)
tmp = np.array([len(hfsyn[kevnm].attrs['comp']) for kevnm in evns])
mcomp = max(tmp)
test = np.array([hfsyn[kevnm].attrs['test'] for kevnm in evns])
if MYID == 0:
  print('current flag and events: ', flag, evns)
if flag == 1: #adj src
  # save to adj.h5 for plot
  if MYID == 0:
    copyfile(fsyn, fadj)
  comm.Barrier()
  hfadj = h5py.File(fadj, 'r+', driver='mpio', comm=MPI.COMM_WORLD)

if flag == 1 or flag ==2: # calculate
  # opt
  misfits_all = np.zeros([nsrc, mcomp])
  energys_all = np.zeros([nsrc, mcomp])
  misfits_part = np.zeros([nsrc, mcomp])
  energys_part = np.zeros([nsrc, mcomp])
  for i in range(MYID, nsrc, NP):
    kevnm = evns[i]
    grpobs = hfobs[kevnm]
    grpsyn = hfsyn[kevnm]
    ntr   = grpobs.attrs['ntr']
    nt    = grpobs.attrs['nt']
    b     = grpobs.attrs['b']
    delta = grpobs.attrs['delta']
    comp  = grpobs.attrs['comp']
    ncomp = len(comp)
    Rx    = grpobs['R'][:,0]
    Rz    = grpobs['R'][:,2]
    Sx    = grpobs['S'][0]
    Sz    = grpobs['S'][2]
    if flag == 1: #adj src
      grpadj = hfadj[kevnm]
    for j in range(ncomp):
      cmp = comp[j]
      obs = grpobs[cmp][()]
      syn = grpsyn[cmp][()]
      ####################################################
      ####################################################
      for itr in range(ntr):
        if Lnorm == 2:
          misfit, energy, adjsrc = misf_wfm(syn[itr,:], obs[itr,:], delta)
          misfits_part[i, j] += misfit
          energys_part[i, j] += energy
        elif Lnorm ==3:
          misfit, adjsrc = misf_tt(syn[itr,:], obs[itr,:], delta)
          misfits_part[i, j] += misfit
        else:
          raise ValueError('implement misfit function later...')
        if flag == 1: #adj src
          grpadj[cmp][itr,:] = adjsrc
      # save adjoint source
      if flag ==1:
        if cmp == 'Vx':
          fout = os.path.join(fdio, kevnm, 'adjx_shot1.su')
        elif cmp == 'Vz':
          fout = os.path.join(fdio, kevnm, 'adjy_shot1.su')
        elif cmp == 'P':
          fout = os.path.join(fdio, kevnm, 'adj_shot1.su')
        adj = grpadj[cmp][()] # forward time
        Rcoord = np.vstack([Rx, Rz]).T
        write_su(fout, adj, [Sx, Sz], Rcoord, delta)
      ####################################################
      ####################################################
  hfobs.close()
  hfsyn.close()
  if flag == 1: #adj src
    hfadj.close()
  # sum up misfit
  comm.Barrier()
  comm.Reduce([misfits_part, MPI.FLOAT], [misfits_all, MPI.FLOAT], op=MPI.SUM, root=0)
  comm.Reduce([energys_part, MPI.FLOAT], [energys_all, MPI.FLOAT], op=MPI.SUM, root=0)
  if MYID == 0:
    # for optim
    res = misfits_all.sum()
    fp = open(fcost,'w')
    fp.write('%.6e\n'% res)
    fp.close()
    # store the result for line search at step = 0
    if flag == 1: #adj src
      misfits_all2 = misfits_all[test==1]
      res = misfits_all2.sum()
      fp = open(fcost+'_0','w')
      fp.write('%.6e\n'% res)
      fp.close()
      # df = pd.DataFrame(data = {'station':evns})
      # df['misfit'] = np.sum(misfits_all,axis=1)
      # df['energy'] = np.sum(energys_all,axis=1)
      # df.to_csv(fmisfit,float_format="%.6e",index=False)
elif flag == 3: # collect
  if MYID == 0:
    print('collect misfit for step = 0...')
    copyfile(fcost+'_0', fcost)
    # df = pd.read_csv(fevn_info)
    # evns = list(df['station'])
    # test = df['test'].to_numpy()
    # df2 = pd.read_csv(fmisfit) # all misfits during cal. adj src
    # evns2 = list(df2['station'])
    # if evns != evns2:
    #   raise ValueError('There are some wrong in function', __file__)
    # misfits = np.array(df2['misfit'][test==1])
    # # for test
    # res = misfits.sum()
    # fp = open(fcost,'w')
    # fp.write('%.6e\n'% res)
    # fp.close()
else:
  raise ValueError('!!! There are some wrong in function !!!', __file__)

print("Finished...", __file__)