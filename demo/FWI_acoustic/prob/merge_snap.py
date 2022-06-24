#!/usr/bin/env python3
"""
Deconvolve the source wavelet [fin to fGv]
Usage:
    python *.py
"""
import os
import sys
import numpy as np

def calc_grid_for_inv(nxg, nyg, nprocx, nprocy, idx,idy):
    nx,ny =nxg//nprocx, nyg//nprocy
    a=np.ones(nx,np.int32)
    for i in range(2, nprocx+1, 1):
        a=np.hstack([a,i*np.ones(nx,np.int32)])
    
    b=a[::idx]
    sumx=np.zeros(nprocx)
    for i in range(1,nprocx+1,1):
        sumx[i-1]=sum(b==i)
    
    a=np.ones(ny,np.int32)
    for i in range(2, nprocy+1, 1):
        a=np.hstack([a,i*np.ones(ny,np.int32)])
    
    b=a[::idy]
    sumy=np.zeros(nprocy)
    for i in range(1,nprocy+1,1):
        sumy[i-1]=sum(b==i)
    
    return sumx, sumy

#Input paras

fjson       =sys.argv[1]
fin         =sys.argv[2]#"snap.bin.p.shot1"
fout        =sys.argv[3]
flag        =int(sys.argv[4])#>0:forward wfd, <0: backward wfd
paras = eval(open(fjson).read())

nxg, nyg=int(paras['NX']),int(paras['NY'])
nprocx, nprocy=int(paras['NPROCX']),int(paras['NPROCY'])
idx,idy=int(paras['IDX']),int(paras['IDY'])
sumx,sumy=calc_grid_for_inv(nxg, nyg, nprocx, nprocy, idx,idy)
fname=fin+".0.0"
size=os.path.getsize(fname)
nxi,nyi=int(sumx[0]),int(sumy[0])
nt=int(size//4/nxi/nyi)

blk_T  = []
for pos1 in range(nprocx):
    Y_T=[]
    for pos2 in range(nprocy):
        fname=fin+".%d.%d"%(pos1,pos2)
        nxi,nyi=int(sumx[pos1]),int(sumy[pos2])
        print("fname, nt,nxi,nyi: ",fname, nt,nxi,nyi)
        if int(os.path.getsize(fname)/nt/nxi/nyi)!=4:
            raise ValueError("check the file size for "+fname)
        wfd=np.fromfile(fname, np.float32).reshape([nt,nxi,nyi])
        if flag<0:
            wfd=wfd[::-1,:,:]
        Y_T.append(wfd)
    if len(Y_T)>0:
        blk_T.append(Y_T)

if len(blk_T)>0:
    np.block(blk_T).tofile(fout)#axis (t, x, z)

print("\nxmovie n1=%d n2=%d < %s loop=1 label1=Y label2=X title=%s\n"% \
     (nyg, nxg,fout,"%g"))

print("Finished...", __file__)
