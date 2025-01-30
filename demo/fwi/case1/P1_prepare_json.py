#!/usr/bin/env python3
"""
Usage:
"""
import os
import sys
import json
import shutil
import numpy as np
import pandas as pd

#Input paras
fdir                    =sys.argv[1]
foptim                  =os.path.join(fdir,'optim.json')
optim                   =eval(open(foptim).read())
fdio                    =optim['fdio']
fevn_info               =optim['fevn_info']
fconf                   =optim['fconf']
nproc                   =optim['nproc']
nxzfd                   =optim['nxzfd']
dxzfd                   =optim['dxzfd']
nt                      =optim['nt']
delta                   =optim['delta']
tshift_adj              =optim['tshift_adj']
freesurface             =optim['freesurface']
pml                     =optim['pml']
pmlvel                  =optim['pmlvel']
comp                    =optim['comp']
epsilon                 =optim['epsilon']

# read station info.
# station,x,z
df = pd.read_csv(fevn_info)
evns = list(df['station'])

num = len(evns)
os.makedirs(os.path.join(fdio, 'comm'), exist_ok=True)
for evn in evns:
    os.makedirs(os.path.join(fdio, evn), exist_ok=True)
    
# write json file into comm 
##################################################
fjson = os.path.join(fconf, 'FD.json')
pars = eval(open(fjson).read())
# temporary configure for FD
#Domain Decomposition
pars['NPROCX'] = str(nproc[0])
pars['NPROCY'] = str(nproc[1])
#space
pars['NX'] = str(nxzfd[0])
pars['NY'] = str(nxzfd[1])
pars['DH'] = str(dxzfd[0])
#boundary
pars['FREE_SURF'] = str(freesurface)
pars['FW'] = str(pml)
pars['VPPML'] = str(pmlvel)
#time
pars['TIME'] = str(nt * delta)
pars['DT'] = str(delta)
#Receiver
if comp == ['P']:
  pars['SEISMO'] = str(2) #pressure
else:
  pars['SEISMO'] = str(1) # Vx, Vy

#gradient
pars['TSHIFT'] = str(tshift_adj)
pars['EPSILON'] = str(epsilon)

# preparation for FD simulation
fp = open(os.path.join(fdio, 'comm', 'test.json'),'w')
fp.write(json.dumps(pars,indent=4))
fp.close()

#adjoint modeling
pars['SEISMO2'] = str(1)
fp = open(os.path.join(fdio, 'comm', 'grad.json'),'w')
fp.write(json.dumps(pars,indent=4))
fp.close()

# copy files
for i in range(num):
    evn = evns[i]
    try:
      shutil.copy(os.path.join(fdio, 'comm', 'test.json'), os.path.join(fdio, evn, 'test.json'))
      shutil.copy(os.path.join(fdio, 'comm', 'grad.json'), os.path.join(fdio, evn, 'grad.json'))
    except IOError as e:
      print("Unable to copy file. %s" % e)
    except:
      print("Unexpected error:", sys.exc_info())

print("Finished...", __file__)