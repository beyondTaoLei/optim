#!/usr/bin/env python3
"""
Usage:
"""
import os
import sys
import obspy
import numpy as np
import pandas as pd
from shutil import copyfile
import matplotlib.pylab as plt
import obspy.signal.filter as flt

#Input paras
fdir          =sys.argv[1]
flag          =int(sys.argv[2]) # 0: copy wavelet, 1: filter wavelet
foptim        =os.path.join(fdir,'optim.json')
optim         =eval(open(foptim).read())
fdio          =optim['fdio']
fsource       =optim['fsource']
fevn_info     =optim['fevn_info']
fc            =optim['fc']#4
order         =optim['order']#4
zero_phase    =optim['zero_phase']

if zero_phase==0:
    zero_phase = False  
else:
    zero_phase = True

# read event info.
df = pd.read_csv(fevn_info)
events = list(df['station'])
nevn = len(events)

################ export source ################
print('\n you will write this wavelet for FD operator \n')
for i in range(nevn):
  evn = events[i]
  fin = os.path.join(fsource, evn+'.su')
  fou = os.path.join(fdio, evn, 'wavelet.su')
  if not os.path.exists(fin):
    raise ValueError("no file: "+fin)
  if flag == 0:
    print("copy file: from %30s to %30s"%(fin, fou))
    copyfile(fin,  fou)
  else:
    print("Butterworth filtering(fc order= %.2f %d ) at (%s --> %s)"%(fc, order, fin, fou))
    seis = obspy.read(fin)
    delta = seis[0].stats.delta
    ntr = seis.count()
    # lowpass filter
    for itr in range(ntr):
      seis[itr].data[:] = flt.lowpass(seis[itr].data, fc, 1.0/delta, order, zerophase=zero_phase)
    seis.write(fou)

# old
trace = obspy.read(fin)[0]
old = trace.data
tr_header = trace.stats.su.trace_header
nt = tr_header.number_of_samples_in_this_trace
delta = tr_header.sample_interval_in_ms_for_this_trace/1.0e6
time = np.arange(nt) * delta
# new
trace = obspy.read(fou)[0]
new = trace.data
print('the output wavelet is: ', fou)
plt.plot(time, old/np.max(old), 'r')
plt.plot(time, new/np.max(new), 'b')
plt.xlabel('time (s)')
plt.ylabel('amplitude')
plt.show()

print("Finished...", __file__)