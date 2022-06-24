#!/usr/bin/env python3
"""
generate initial model and input json file for Rosenbrock problem
Usage:
    python S0_prepare.py optim_PSTD
    python S0_prepare.py optim_PNLCG
    python S0_prepare.py optim_LBFGS
    python S0_prepare.py optim_TRN
"""

import os
import sys
import numpy as np
import pandas as pd
from shutil import copyfile

def split_jobs(jobs, cell, ncores):
    """
    split jobs into subjobs for one lsf file
    """
    nsrc=len(jobs)
    num_in_group = ncores//cell
    subjobs=[jobs[i:i+num_in_group] for i in range(0, nsrc, num_in_group)]
    return subjobs

def generate_lsf(jobnm, joblist, queuenm, ncores, wtime, headers, parallel=1):
    """
    submit jobs with lsf system
    """
    fp = open(jobnm,'w')
    a=fp.write("#!/bin/bash\n")
    a=fp.write("#BSUB -q %s\n"%queuenm)
    a=fp.write("#BSUB -n %d\n"%ncores)
    a=fp.write("#BSUB -W %s\n"%wtime)
    a=fp.write("#BSUB -e %J.err\n")
    a=fp.write("#BSUB -o %J.out\n")
    #a=fp.write('#BSUB -R "span[ptile=40]"\n\n')
    a=fp.write("\n")
    a=fp.write("module load mpi/intel/2018.4\n")
    a=fp.write("module load intel/2018.4\n")
    a=fp.write("module load fftw/3.3.8\n")
    a=fp.write("module load matlab/2017b\n\n")
    for s in headers:
        a=fp.write("%s\n"%s)
    a=fp.write("\n")
    
    if parallel == 1:
        for s in joblist[:-1]:
            a=fp.write('%s&\n'% s) 
        a=fp.write('%s\n'% joblist[-1])
    else:
        for s in joblist:
            a=fp.write('%s\n'% s)
    fp.close()

#Input paras
mpiexec     ='/share/intel/2018u4/compilers_and_libraries_2018.5.274/linux/mpi/intel64/bin/mpiexec'
prjroot     ='/data/ess-leit/Marmousi2_new'
spfroot     ='/work/ess-leit/seisplatform'
FDEXE       ='/work/ess-leit/seisplatform/bin/sofi2Dac'
ncores      =40
NP          =4
TESTSHOTINC =3
NT          =6500
DT          =1.0e-03

headers =[
    "prjroot='%s'"%prjroot,
    "FDEXE='%s'"%FDEXE,
    "mpiexec='%s'"%mpiexec,
    "NP=%d"%NP
    ]
queue_we = 'short'
wtime = '2:00'
dirjob= 'job'
#input files
fsrc    ='list/source.dat'
frec    ='list/receiver.dat'
fsrcname='list/eventname.list'
frecname='list/stationname.list'
ffdjson ='list/FW.json'
ffdjson2='list/GRAD.json'

#intermediate files
fevn_info  ='list/eventinfo.csv'
fstn_info  ='list/stationinfo.csv'
fdio       ='fdio'
fsource    ='source'

# make directory
os.makedirs(fdio, exist_ok=True)
os.makedirs(fsource, exist_ok=True)

#wavelet
cmd0=os.path.join(spfroot,'bin','write_wavelet')\
    +' nt='+str(NT)+' delta='+str(DT)\
    +' fsrc='+fsrc+' fpath=source'
print(cmd0)
#p=os.system(cmd0)
#if p!=0:
    #raise ValueError('you fail to run write_wavelet')

#source 
#coordinate
Scoord=np.loadtxt(fsrc, skiprows=1)
#event name
a=pd.read_csv(fsrcname,header=None, sep=r"[| ^]", engine='python',names=['id','stn'])
shot_list=list(a['stn'])
d={'station':shot_list}
df=pd.DataFrame(data=d)
df["X"]     =Scoord[:,0]
df["tmp"]   =Scoord[:,1]
df["Y"]     =Scoord[:,2]
df["tshift"]=Scoord[:,3]
df["fc"]    =Scoord[:,4]
df["amp"]   =Scoord[:,5]
testflag=np.zeros(len(shot_list),np.int32)
testflag[::TESTSHOTINC]=1
df["test"]=testflag
del df["tmp"]
df.to_csv(fevn_info,index=False)

#receiver
#coordinate 2
Scoord=np.loadtxt(frec)
#station name
a=pd.read_csv(frecname,header=None, sep=r"[| ^]", engine='python',names=['id','stn'])
recv_list=list(a['stn'])
d={'station':recv_list}
df=pd.DataFrame(data=d)
df["X"] =Scoord[:,0]
df["Y"] =Scoord[:,1]
df.to_csv(fstn_info,index=False)

#FD work directory
#copy model receiver.dat input.json into each work directory
print(70*'*')
print('copy [source.dat|receiver.dat|FW.json|GRAD.json] into fdio/comm \n')
print(70*'*')
os.makedirs(os.path.join(fdio, 'comm'), exist_ok=True)
copyfile(fsrc,      os.path.join(fdio, 'comm','source.dat'))
copyfile(frec,      os.path.join(fdio, 'comm','receiver.dat'))
copyfile(ffdjson,   os.path.join(fdio, 'comm','FW.json'))
copyfile(ffdjson2,  os.path.join(fdio, 'comm','GRAD.json'))

for i in range(len(shot_list)):
    os.makedirs(os.path.join(fdio, shot_list[i]), exist_ok=True)
    fin=os.path.join(fsource, 'wavelet_shot'+str(i+1)+'.su')
    copyfile(fin,  os.path.join(fdio, shot_list[i], 'wavelet.su'))
    
#generate all commands
paths=['cd '+os.path.join('${prjroot}', 'fdio', evt) for evt in shot_list]
function_forward = '${mpiexec}' + ' -np ${NP} ${FDEXE}'+ ' ../comm/FW.json'
function_data_process = 'mv seis_p_shot1.su obs.su'
function_modelling_shot = [function_forward, function_data_process]

jobs=[]
for isrc in range(len(shot_list)):
    cmd=[paths[isrc]] + function_modelling_shot
    jobs.append('&& '.join(cmd))
subjobs = split_jobs(jobs, NP, ncores)
for i in range(len(subjobs)):
    jobnm = os.path.join(dirjob, 'S00_part%d.lsf'%(i+1))
    generate_lsf(jobnm, subjobs[i], queue_we, ncores, wtime, headers, 1)

print('*'*60)
print("Finished...", __file__)
#####################