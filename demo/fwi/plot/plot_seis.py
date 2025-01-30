#!/usr/bin/env python3
import sys
import h5py
import numpy as np
import matplotlib.pylab as plt
import matplotlib
# matplotlib.use('TKAgg') # for Qt error

print('plot_seis.py data.h5 evnm cmp vmax')
fin   =sys.argv[1] # cc.h5
evnm  =sys.argv[2] # sz001
cmp   =sys.argv[3] # Vz
if len(sys.argv)==5:
  vmax  =float(sys.argv[4])
  vmin  =-vmax

# figure name to save
fignm = fin[:-3]+'_'+evnm+'.png'
# figure size to save
figsize = [3,5]
figsize = [4,2]
# figure resolution to save
figdpi  = 150
##########################################
hf = h5py.File(fin, 'r')
grp = hf[evnm]
data = grp[cmp][()]
ntr, nt = data.shape
b  = grp.attrs['b']
# nt  = grp.attrs['nt']
delta = grp.attrs['delta']
time = b + np.arange(nt) * delta
hf.close()
if ntr == 1:
  data = np.tile(data, (5,1))
  ntr = 5
print('name, cmp, ntr, nt, data.max()',evnm, cmp, ntr, nt, data.max())
if len(sys.argv)==4:
  abs1 = abs(data)
  vmax = np.percentile(abs1, 90)
  vmin  =-vmax

# show
fig, ax = plt.subplots(dpi=figdpi,figsize=(figsize[0],figsize[1])) #jet, gray_r
plt.imshow(data.T, cmap='gray_r', aspect='auto',
          extent=[1, ntr, time[-1], time[0]], vmin=vmin, vmax=vmax)

plt.colorbar()
# plt.plot(np.arange(ntr), np.zeros(ntr),'r')
ax.set(title=fin[:-3]+':'+evnm, xlabel='Trace No.', ylabel='Time (s)')
plt.savefig(fignm)
plt.show()

# plt.plot(time, data[1,:])
# plt.show()

# print("Finished...", __file__)