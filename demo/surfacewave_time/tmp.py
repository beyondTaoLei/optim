import os
import sys
import json

fdir        =sys.argv[1]
iter0       =int(sys.argv[2])
foptim      =os.path.join(fdir,'optim.json')
optim       =eval(open(foptim).read())
optim['iter0'] = iter0

fout = open(foptim,'w')
fout.write(json.dumps(optim,indent=4))
fout.close()