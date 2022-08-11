#!/usr/bin/env python3
import os
import sys
import shutil

#start here    
#Input paras
fdir        =sys.argv[1] #optim
idx1        =sys.argv[2]#'6'
foptim      =os.path.join(fdir,'optim.json')
optim       =eval(open(foptim).read())
iterc       =optim['iterc']

fin = os.path.join(fdir, 'iter'+str(iterc), 'L_curve_'+str(idx1)+'.bin')
fout = os.path.join(fdir, 'iter'+str(iterc), 'desc.bin')
print("copy "+fin+ " to "+fout)
shutil.copy(fin, fout)

flog = os.path.join(fdir, 'iter'+str(iterc), 'L_curve.log')
fp = open(flog,'w')
fp.write("copy "+fin+ " to "+fout+"\n")
fp.close()

print("Finished...", __file__)