#!/usr/bin/env python3
"""
Usage:
        
"""

import os
import sys
import json

#Input paras
fdir        =sys.argv[1]
key         =sys.argv[2]
if len(sys.argv) == 4:
    value_new   =sys.argv[3]
foptim      =os.path.join(fdir,'optim.json')
optim       =eval(open(foptim).read())
value_old   =optim[key]
type1       =type(optim[key])

print('*'*30)
print('key : ', key)
print('type: ', type1)
print('old : ', value_old)
if len(sys.argv) == 4:
    optim[key]  =type1(value_new)
    print('new : ', value_new)
print('*'*30)

# write work infomation
if len(sys.argv) == 4:
    foptim = os.path.join(fdir,'optim.json')
    fout = open(foptim,'w')
    fout.write(json.dumps(optim,indent=4))
    fout.close()
    print('*'*60+'\n*Output file: '+foptim)

print('Finished...', __file__)