#!/usr/bin/env python3
"""
useful functions
Usage:
    [from util import fun] in another *.py file
https://stackoverflow.com/questions/20309456/call-a-function-from-another-file
https://docs.obspy.org/master/packages/autogen/obspy.io.segy.header.TRACE_HEADER_FORMAT.html#obspy.io.segy.header.TRACE_HEADER_FORMAT
"""
import os
import math
import obspy
import shutil
import numpy as np
import pandas as pd
from scipy.signal.windows import tukey
from scipy.ndimage import gaussian_filter
from obspy.io.segy.segy import SEGYTraceHeader, SEGYBinaryFileHeader

class Empty:
    def __init__(self):
      self.aa = 0

def write_su(fnm, data, Sxy, Rcoord, delta):
  """
  write data to SU file
  """
  ntr, nt = data.shape
  # generate the wavelets, and save to segy-format file
  scalel = 1000.0
  scalco = 1000.0
  # create Stream
  stream = obspy.Stream()
  for itr in range(ntr):
    trace = obspy.Trace(data=data[itr, :])
    # Attributes in trace.stats will overwrite everything in
    # trace.stats.segy.trace_header
    trace.stats.starttime = obspy.UTCDateTime(1970,1,1,00,00,00)
    trace.stats.delta = delta
    # format: obspy.io.segy.header.TRACE_HEADER_FORMAT
    trace.stats.su = {}
    trace.stats.su.trace_header = SEGYTraceHeader()
    tr_header = trace.stats.su.trace_header
    tr_header.trace_sequence_number_within_line = itr
    tr_header.trace_identification_code = 1 # trace identification code: 1=seismic
    tr_header.receiver_group_elevation = int(Sxy[1] * scalel) # not use
    tr_header.water_depth_at_source = 1000 # not use
    tr_header.scalar_to_be_applied_to_all_elevations_and_depths = - int(scalel)
    tr_header.scalar_to_be_applied_to_all_coordinates = -int(scalco)
    tr_header.source_coordinate_x   = int(Sxy[0] * scalco) #Sx
    tr_header.source_coordinate_y   = int(Sxy[1] * scalco) #Sy
    tr_header.group_coordinate_x    = int(Rcoord[itr, 0] * scalco) #Rx=Sx
    tr_header.group_coordinate_y    = int(Rcoord[itr, 1] * scalco) #Ry=Sy
    tr_header.number_of_samples_in_this_trace = nt 
    tr_header.sample_interval_in_ms_for_this_trace = int(delta * 1.0e6) # sample interval in micro-seconds
    # Add trace to stream
    aa = stream.append(trace)
  # header
  stream.stats = obspy.core.AttribDict()
  stream.stats.textual_file_header = 'Textual Header!'
  stream.stats.binary_file_header = SEGYBinaryFileHeader()
  stream.stats.binary_file_header.trace_sorting_code = 5
  stream.write(fnm, format='SU', byteorder='<')
  #print('suswapbytes <*.su format=0 ns=6667 |supsimage perc=0.96 >1.eps')

def taper2d(n2, n1, iscoord, r1, r2):
  """
  taper source
  """
  window1d = 1.0 - np.abs(tukey(2*r2+1, (r2-r1)/r2))
  x = np.arange(2*r2+1) - r2
  X = np.tile(x, [len(x),1]).T
  dist = np.sqrt(X**2 + X.T**2)
  window2d = np.interp(dist, x, window1d)
  window2d = np.where(dist<=r2, window2d, 1.0)
  nw2, nw1 = window2d.shape
  taper = np.ones([n2+2*r2, n1+2*r2],np.float32)
  for i, j in iscoord:
    taper[i:i+nw2,j:j+nw1] *= window2d
  #kernel part
  return taper[r2:r2+n2,r2:r2+n1]

def taper_shot(x, z, srcpos, r1, r2):
  # find index
  iarr = np.zeros_like(srcpos, np.int32)
  num, tmp = srcpos.shape
  for i in range(num):
    x0, z0 = srcpos[i]
    dist = (x0 - x)**2 + (z0 - z)**2
    iarr[i,:] = np.unravel_index(dist.argmin(), dist.shape) #first iz, then ix
  # apply taper
  nz, nx = x.shape
  taperzx = taper2d(nz, nx, iarr, r1, r2)
  return taperzx

def ttaper(t, t1, t2, nalpha=20):
  """
  select one phase
  """
  nt = len(t)
  taper = np.zeros_like(t, dtype = np.float32)
  it1 = np.argmin(abs(t-t1))
  it2 = np.argmin(abs(t-t2))
  nt2 = it2 - it1 + 1
  taper[it1:it2+1] = tukey(nt2, nalpha/nt)
  return taper

def ttaperL(t, t1, nalpha=20):
  """
  select one phase
  """
  nt = len(t)
  taper = np.zeros_like(t, dtype = np.float32)
  it1 = np.argmin(abs(t-t1))
  a = tukey(2*nalpha, 1.0)
  taper[it1:it1+nalpha] = a[:nalpha]
  taper[it1+nalpha:] = 1.0
  return taper

def ttaperR(t, t1, nalpha=20):
  """
  select one phase
  """
  nt = len(t)
  taper = np.zeros_like(t, dtype = np.float32)
  it1 = np.argmin(abs(t-t1))
  a = tukey(2*nalpha, 1.0)
  taper[it1-nalpha:it1] = a[nalpha:]
  taper[:it1-nalpha] = 1.0
  return taper

# Function to find the smallest power of 2
# greater than or equal to n
def nextpow2(N):
    # Calculate log2 of N
    a = int(math.log2(N))
    # If 2^a is equal to N, return N
    if 2**a == N:
        return N
    # Return 2^(a + 1)
    return 2**(a + 1)

def fk_filter(nx, nt, dx, dt, cmin, cmax, sigma=3.0):
  np.seterr(divide='ignore', invalid='ignore')
  nk = 4*(2^nextpow2(nx))
  nf = 4*(2^nextpow2(nt))
  # freq = np.arange(-nf//2+1, nf//2+1)/nf/dt
  # kx = np.arange(-nk//2+1, nk//2+1)/nk/dx
  freq = np.fft.fftshift(np.fft.fftfreq(nf, d=dt))
  kx = np.fft.fftshift(np.fft.fftfreq(nk, d=dx))
  kxs = np.tile(kx,[nf,1]).T
  freqs = np.tile(freq,[nk,1])
  v = np.where(abs(kxs) > 0.1*(kx[1]-kx[0]), abs(freqs/kxs), 10000.0)
  filt = np.where((v>=cmin)&(v<=cmax), 0.0, 1.0)
  # smooth fk-filter using 2D Hamming window
  filt[:] = gaussian_filter(filt, sigma=sigma)
  return kx, freq, filt

def linefitting(Xs, Ys):
  '''
  Least square method line fitting Y = a*X + b
  '''
  N = len(Xs)
  l1 = np.sum(Xs)
  l2 = np.sum(Ys)
  l3 = np.sum(Xs*Ys)
  l4 = np.sum(Xs*Xs)
  a = (N * l3 - l1 * l2) / (N * l4 - l1 * l1)
  b = (l4* l2 - l1 * l3) / (N * l4 - l1 * l1)
  ave = l2 / N
  return a, b, ave

def fun(a,b):
    return a+b
