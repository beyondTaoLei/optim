#!/usr/bin/env python3
import os
import sys
import numpy as np
import pandas as pd
from util import FSM_rays
from util import fastsweep
from util import extract_record
from util import sourceInitiation

#input paras
fdir        =sys.argv[1]
foptim      =os.path.join(fdir,'optim.json')
optim       =eval(open(foptim).read())
ox          =optim['o2g']
oy          =optim['o1g']
nx          =optim['n2g']
ny          =optim['n1g']
dx          =optim['d2g']
dy          =optim['d1g']
fmodfd      =optim['fmodfd']
fmod        =os.path.join('../..', fmodfd)
frec        ='receiver.csv'
ftime       ='time.csv'
fraypath    ='raypath.dat'
ftimetable  ='timefield.bin'
itermaxFSM  =3
ixy         =5
eps         =1.0e-08

#interpretation
xvel=np.arange(nx)*dx+ox
yvel=np.arange(ny)*dy+oy
if not os.path.exists(frec):
    print('No receiver file at '+ os.getcwd())
    exit(1)
df=pd.read_csv(frec)
coordR=np.vstack([df['stx'], df['sty']]).transpose()
#print(coordR)
if 'sxy' in optim.keys():
    sxy =optim['sxy']#"2500.0,2500.0"
    sx, sy=[float(i) for i in sxy.split(',')]
else:
    sx, sy=df['evx'][0],df['evy'][0]
valinfty    =1000.0
slow=np.fromfile(fmod, np.float32).reshape([nx,ny]).astype(np.float64)
slow[:]=1.0/slow
#############################
# initialize the time field
#############################
timetable = valinfty*np.ones([nx, ny], np.float64)
source = np.zeros([nx, ny], np.float32)
sourceInitiation(sx, sy, xvel, yvel,timetable, source, slow, ixy)
iteration = 0
gap = 2.0*eps
while iteration < itermaxFSM and gap > eps:
    gap=fastsweep(0, nx-1, 0, ny-1, nx, ny, dx, dy, slow, timetable, source)
    gap=fastsweep(nx-1, 0, 0, ny-1, nx, ny, dx, dy, slow, timetable, source)
    gap=fastsweep(0, nx-1, ny-1, 0, nx, ny, dx, dy, slow, timetable, source)
    gap=fastsweep(nx-1, 0, ny-1, 0, nx, ny, dx, dy, slow, timetable, source)
    iteration = iteration + 1
print('Iteration ', iteration, ' and gap is ',gap)

timetable.astype(np.float32).tofile(ftimetable)
time=extract_record(coordR, xvel, yvel,timetable)
df['time']=time
df.to_csv(ftime, float_format="%.4f", index=False)

# ray path
Tx=np.diff(timetable, axis=0) #out[i] = a[i+1] - a[i] 
Tx=np.vstack([Tx, Tx[-1,:]])/dx
Ty=np.diff(timetable, axis=1)
Ty=np.hstack([Ty, np.atleast_2d(Ty[:,-1]).T])/dy
FSM_rays(fraypath, xvel, yvel, sx, sy, coordR, 1.0/slow[:], timetable, Tx, Ty, dx)

#print("Finished...", __file__)    
    
        
