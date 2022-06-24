#!/usr/bin/env python3
"""
copy extern files to optim file system at each iteration
"""
import os, sys, json
from shutil import copyfile

#Input paras
fdir    =sys.argv[1] #'optim'
iterc   =int(sys.argv[2])
foptim  =os.path.join(fdir,'optim.json')
optim   =eval(open(foptim).read())
fmod    =optim['fmod']
fcost   =optim['fcost']
fgrad   =optim['fgrad']
iter0   =optim['iter0']
fiter   ='iter'+str(iterc)

#save arrays to optim file system
os.makedirs(os.path.join(fdir, fiter),exist_ok=True)
copyfile(fmod,  os.path.join(fdir, fiter, 'mod.bin'))
copyfile(fcost, os.path.join(fdir, fiter, 'fcost.dat'))
copyfile(fgrad, os.path.join(fdir, fiter, 'grad.bin'))

#update work infomation
optim['iterc'] =iterc
if iter0==iterc: #special for l-BFGS
    optim['cpt_lbfgs'] =0
fout =open(foptim,'w')
fout.write(json.dumps(optim,indent=4))
fout.close()

#print("Finished...", __file__)
  
   
