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
import obspy.signal.filter as flt
from scipy.ndimage import gaussian_filter
from misfit_funs import misf_wfm, misf_tt
from util import ttaper, ttaperL, fk_filter

#initialize MPI
comm=MPI.COMM_WORLD
NP = comm.Get_size()
MYID = comm.Get_rank()

#Input paras
fdir              =sys.argv[1]
foptim            =os.path.join(fdir,'optim.json')
optim             =eval(open(foptim).read())
fdio              =optim['fdio']
fdata             =optim['fdata']
dirmute           =optim['dirmute']
dirmute_taperlen  =optim['dirmute_taperlen']
dirmute_t0        =optim['dirmute_t0']
dirmute_lower     =optim['dirmute_lower']
dirmute_upper     =optim['dirmute_upper']
timefilt          =optim['timefilt']
fc                =optim['fc']
order             =optim['order']
zero_phase        =optim['zero_phase']
timewin           =optim['timewin']
timewin_lower     =optim['timewin_lower']
timewin_upper     =optim['timewin_upper']
trkill            =optim['trkill']
trkill_lower      =optim['trkill_lower']
trkill_upper      =optim['trkill_upper']
velocity          =optim['velocity']

#parse inputs
if zero_phase==0:
  zero_phase = False  
else:
  zero_phase = True

#parse inputs
fin = os.path.join(fdata, 'syn0.h5')
fou = os.path.join(fdata, 'syn.h5')
# copy
if MYID == 0:
  copyfile(fin, fou)
comm.Barrier()
hfin = h5py.File(fin, 'r', driver='mpio', comm=MPI.COMM_WORLD)
hfou = h5py.File(fou, 'r+', driver='mpio', comm=MPI.COMM_WORLD)
if dirmute == 3:
  fmute = os.path.join(fdata, 'mute.h5')
  hfmt = h5py.File(fmute, 'r', driver='mpio', comm=MPI.COMM_WORLD)
evns = list(hfin.keys())
nsrc = len(evns)
if MYID == 0:
  print('current events: ', evns)
comm.Barrier()

# opt
ntr0, nt0 = 0, 0
for i in range(MYID, nsrc, NP):
  kevnm = evns[i]
  grpin = hfin[kevnm]
  ntr   = grpin.attrs['ntr']
  nt    = grpin.attrs['nt']
  b     = grpin.attrs['b']
  delta = grpin.attrs['delta']
  comp  = grpin.attrs['comp']
  Sx    = grpin['S'][0]
  Sz    = grpin['S'][2]
  Rx    = grpin['R'][:,0]
  Rz    = grpin['R'][:,2]
  time = b + np.arange(nt) * delta
  dist = np.sqrt((Sx-Rx)**2+(Sz-Rz)**2)
  ncomp = len(comp)
  for j in range(ncomp):
    cmp = comp[j]
    print(MYID, kevnm, cmp)
    data  = grpin[cmp][()]
    # integration of data
    if velocity == 0:
      for itr in range(ntr):
        data[itr,:] = np.cumsum(data[itr,:]) * delta
    # muting direct waves
    if dirmute == 1:
      for itr in range(ntr):
        t0 = dist[itr] / dirmute_lower + dirmute_t0
        nalpha = int(dirmute_taperlen/delta)
        taperd = ttaperL(time, t0, nalpha)
        data[itr,:] *= taperd
    elif dirmute == 2:
      if ntr != ntr0 or nt != nt0:
        dx = abs(Rx[1] - Rx[0])
        kx, freq, filt = fk_filter(ntr, nt, dx, -delta, dirmute_lower, dirmute_upper, 3.0)
        filto = np.fft.ifftshift(filt)
        d2 = np.zeros_like(filt, np.float32)
        ntr0, nt0 = ntr, nt
      d2[:ntr, :nt] = data
      tmp2 = np.fft.ifft2(np.fft.fft2(d2) * filto)
      data[:] = tmp2[:ntr, :nt].real
    elif dirmute == 3:
      grpmt = hfmt[kevnm] 
      tapermt = grpmt['win'][()]
      data *= tapermt
    # filtering trace by trace
    if timefilt == 1:
      for itr in range(ntr):
        data[itr,:] = flt.lowpass(data[itr,:], fc, 1.0/delta, order, zerophase=zero_phase)
    # trace killing
    if trkill == 1:
      # ntr maybe different for each shot
      trcoeff =np.where((dist>=trkill_lower)&(dist<=trkill_upper), 1.0, 0.0) # np.float64
      trcoeff[:] = gaussian_filter(trcoeff, sigma=1.0)
      for itr in range(ntr):
        data[itr,:] *= trcoeff[itr] 
    # applying time window (nt is same for all shots)
    if timewin == 1:
      if i == MYID:
        tapert = ttaper(time, timewin_lower, timewin_upper, 30)
        taper2d = np.tile(tapert,[ntr,1])
      data *= taper2d
    # normalization
    # envelope
    # rewrite
    grpou = hfou[kevnm]
    grpou[cmp][:] = data

hfin.close()
hfou.close()
if dirmute == 3:
  hfmt.close()

print("Finished...", __file__)
