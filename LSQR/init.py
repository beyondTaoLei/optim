#!/usr/bin/env python3
"""
copy extern files to optim file system at each iteration
"""
import os
import sys
import json
from shutil import copyfile

#Input paras
fdir    =sys.argv[1] #optim
iterc   =int(sys.argv[2])
foptim  =os.path.join(fdir,'optim.json')
optim   =eval(open(foptim).read())
fmodfd  =optim['fmodfd']
fmod    =optim['fmod']
fmodref =optim['fmodref']
fdiff   =optim['fdiff']
flag    =optim['modeltype']#1:smoothest, 2:flattest, 3:smallest
fiter   ='iter'+str(iterc)
if 'fcost' in optim:
    fcost   =optim['fcost']

#save arrays to optim file system
os.makedirs(os.path.join(fdir, fiter),exist_ok=True)
copyfile(fmodfd,os.path.join(fdir, fiter, 'modg.bin'))
copyfile(fmod,  os.path.join(fdir, fiter, 'mod.bin'))
copyfile(fdiff, os.path.join(fdir, fiter, 'diff.dat'))
if 'fcost' in optim:
    copyfile(fcost, os.path.join(fdir, fiter, 'fcost.dat'))
if flag==3:
    copyfile(fmodref,  os.path.join(fdir, fiter, 'ref.bin'))

#update work infomation
optim['iterc'] =iterc
fout =open(foptim,'w')
fout.write(json.dumps(optim,indent=4))
fout.close()

#print("Finished...", __file__)