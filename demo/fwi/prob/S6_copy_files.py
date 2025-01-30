#/home/tao/seis_software
import os
import sys
import shutil

#Input paras
fdir        =sys.argv[1]
foptim      =os.path.join(fdir,'optim.json')
optim       =eval(open(foptim).read())
iterc       =optim['iterc']
fiter       ='iter'+str(iterc)

shutil.copyfile('data/syn_bk.h5',       os.path.join(fdir, fiter, 'syn.h5'))
shutil.copyfile(foptim,                 os.path.join(fdir, fiter, 'optim.json'))

print("Finished...", __file__)        