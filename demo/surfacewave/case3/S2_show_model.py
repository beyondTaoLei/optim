#!/usr/bin/env python3
import os
import sys
import json
import numpy as np
import matplotlib.pylab as plt

fdir        =sys.argv[1]
iterc       =sys.argv[2]
depth       =float(sys.argv[3])
foptim      =os.path.join(fdir,'optim.json')
optim       =eval(open(foptim).read())
n3          =optim['n3g']
n2          =optim['n2g']
n1          =optim['n1g']
d3          =optim['d3g']
d2          =optim['d2g']
d1          =optim['d1g']
o3          =optim['o3g']
o2          =optim['o2g']
o1          =optim['o1g']
X = np.arange(n3)*d3 + o3
Y = np.arange(n2)*d2 + o2
Z = np.arange(n1)*d1 + o1

fnm = os.path.join(fdir, 'iter'+iterc, 'mod.bin')
mod = np.fromfile(fnm, np.float32).reshape([n3, n2, n1])

iz = int((depth-o1)//d1)
slice1 = mod[:, :, iz]

x1 = X[0]
x2 = X[-1]
y1 = Y[0]
y2 = Y[-1]
vmin = int(1.05*slice1.min())
vmax = int(0.95*slice1.max())
title1 = 'Z' + str(depth) + 'm'

params = {"figure.figsize":(4,7),
          'axes.labelsize': 14,
          'axes.titlesize': 16,
          'xtick.labelsize': 12,
          'ytick.labelsize': 12,
          'legend.fontsize': 14}
plt.rcParams.update(params)
fig, ax=plt.subplots()
im=ax.pcolormesh(
            X,
            Y, 
            slice1.T,
            cmap='jet_r',
            shading='auto',
            vmin=vmin,
            vmax=vmax)
ax.set_xlim([x1, x2])
ax.set_ylim([y1, y2])
ax.set_xlabel('X (km)')
ax.set_ylabel('Y (km)')
ax.tick_params('both')
plt.title(title1)
###############################
ax.grid(True)
fig.colorbar(im)

#df=pd.read_csv('../list/stationinfo.csv')
#dfx=np.array(df['stx'])/1000
#dfy=np.array(df['sty'])/1000
#ax2 = ax.twinx()
#ax2.scatter(dfx, dfy,s=10,c='w', marker="^")
#ax2.set_ylim([y1, y2])
##for i in range(0, len(dfx),1):
    ##ax2.text(dfx[i],dfy[i],"%s"%(df['station'][i]))
#plt.tight_layout()
#plt.title('Depth '+str(depth)+' km')
plt.savefig('fig/'+title1+'.png')
#plt.savefig('../fig/'+title1+'_2.png')
plt.show()




