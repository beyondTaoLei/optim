#!/usr/bin/env python3
"""
generate initial model and input json file for Rosenbrock problem
Usage:

"""

import os
import sys
import json
import numpy as np
import pandas as pd

os.makedirs('list', exist_ok=True)

# source
ntr = 160
Rcoord = np.zeros([ntr, 2], np.float32)
Rcoord[:,0] = np.arange(ntr) * 100 + 500
Rcoord[:,1] = 200.0
stn_list = ['sz%03d'%(i+1) for i in range(ntr)]

Scoord = Rcoord[::4,:]
evn_list = stn_list[::4]
nsrc = len(evn_list)

# srcinfo
df = pd.DataFrame(data={'station':evn_list})
df['x'] = Scoord[:,0]
df['y'] = Scoord[:,1]
df.to_csv('list/src.csv',float_format="%.2f",index=False)

#kevnm,kstnm,evx,evy,stx,sty,dist
d = {'kevnm':[evn_list[0]]*ntr}
df = pd.DataFrame(data=d)
df["kstnm"] = stn_list
df["evx"] = Scoord[0,0]
df["evy"] = Scoord[0,1]
df["stx"] = Rcoord[:,0]
df["sty"] = Rcoord[:,1]

os.makedirs('fdio', exist_ok=True)
for isrc in range(nsrc):
    fsrc0 = os.path.join('fdio', evn_list[isrc])
    os.makedirs(fsrc0, exist_ok=True)
    df["kevnm"] = [evn_list[isrc]]*ntr
    df["evx"] = Scoord[isrc,0]
    df["evy"] = Scoord[isrc,1]
    df.to_csv(os.path.join(fsrc0,'receiver.csv'),float_format="%.2f",index=False)

print('Finished...', __file__)
#####################