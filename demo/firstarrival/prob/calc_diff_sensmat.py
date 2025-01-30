#!/usr/bin/env python3
import os
import sys
import numpy as np
import pandas as pd
from scipy import sparse
from numpy import linalg as LA
from scipy.ndimage import gaussian_filter
from util import raypath2senstivity

#input paras
fdir        =sys.argv[1]
foptim      =os.path.join(fdir,'optim.json')
optim       =eval(open(foptim).read())
nx, ny      =optim['nxyfd']
ox, oy      =optim['oxyfd']
dx, dy      =optim['dxyfd']
incx, incy  =optim['incxy']
fsrc        =optim['fsrc']
fcost       =optim['fcost']
fdiff       =optim['fdiff']
fsensmat    =optim['fsensmat']
ftimeobs    ='obs.csv'
ftime       ='time.csv'
fraypath    ='raypath.dat'

#interpretation
xvel = np.arange(nx)*dx+ox
yvel = np.arange(ny)*dy+oy
nxi = len(np.arange(nx)[::incx])
nyi = len(np.arange(ny)[::incy])

# read source info.
evn_info=pd.read_csv(fsrc)
shot_list=list(evn_info['station'])

# read time
times = {}
for evt in shot_list:
    fobs = os.path.join('fdio', evt, ftimeobs)
    fsyn = os.path.join('fdio', evt, ftime)
    dfobs = pd.read_csv(fobs)
    dfsyn = pd.read_csv(fsyn)
    tobs = np.array(dfobs['time'], np.float32)
    tsyn = np.array(dfsyn['time'], np.float32)
    diff=np.where((tsyn>0.0) &(tobs>0.0), tobs-tsyn, -1000.0)
    dfsyn.rename(columns={'time':'syn'}, inplace=True)
    dfsyn['obs'] = tobs
    dfsyn['diff'] = diff
    times[str(evt)] = dfsyn

# save the time difference
df = pd.concat([times[str(evt)] for evt in shot_list])
df.to_csv(fdiff[:-4]+'.csv', float_format="%.4f", index=False)
b=np.array(df['diff'])

# read raypath
raypaths_list = []
for evt in shot_list:
    fin = os.path.join('fdio', evt, fraypath)
    npoints_arr = np.loadtxt(fin, max_rows=1)
    raypaths = np.loadtxt(fin, skiprows=1)
    sum1 = np.cumsum(npoints_arr)
    idx = np.hstack([[0],sum1]).astype(np.int32)
    for i in range(len(idx)-1):
        raypaths_list.append(raypaths[idx[i]:idx[i+1]])

# filtering
if len(b) != len(raypaths_list):
    raise ValueError("check your forward module!")

idxs = np.where(b>-900.0)[0]
b = b[idxs]
nxg2, nyg2 = nxi*incx, nyi*incy
mat2d = np.zeros([nxg2, nyg2], np.float32)
res = []
for j in idxs:
    #print(j)
    raypath2senstivity(raypaths_list[j], xvel, yvel, mat2d)
    arr2d = mat2d[:nxg2, :nyg2]
    mat2dblk = arr2d.reshape(nxi, incx, nyi, incy).sum(axis=(1, 3))
    mat2dblk = gaussian_filter(mat2dblk, sigma=3.0*incx/np.sqrt(8.0))
    cscmat = sparse.csc_matrix(mat2dblk.flatten())
    res.append(cscmat)

# save sensitivity matrix with format of csc_matrix
sen_mat = sparse.vstack(res).tocsc()
sparse.save_npz(fsensmat, sen_mat)
np.savetxt(fdiff, b, fmt='%.4e')

misfit=0.5*LA.norm(b)**2
fp = open(fcost,'w')
fp.write('%.6e\n'% misfit)
fp.close()

#print("Finished...", __file__)
    
        
