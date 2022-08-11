#!/usr/bin/env python3
"""
calculates descent direction
python True=1 False=0
"""
import os
import sys
import json
import numpy as np
from CG_test import test_singularity, test_negative_curv
from CG_test import test_negative_curv2, test_truncation, record_min_residual
from forcing_term_TRN import forcing_term_TRN
# convert dictionary into a instance of class
#https://www.geeksforgeeks.org/how-to-change-a-dictionary-into-a-class/
class Dict2Class(object):  
    def __init__(self, my_dict):
        for key in my_dict:
            setattr(self, key, my_dict[key])

def load_worktn(fwork, optim):
    if optim.nwork!=10:
        raise ValueError("change nwork firstly")
    worktn=np.fromfile(fwork, np.float32).reshape([optim.nwork, optim.n])
    optim.worktn        =worktn
    optim.mod           =worktn[0,:]
    optim.grad          =worktn[1,:]
    optim.descent       =worktn[2,:]
    optim.descent_prev  =worktn[3,:]
    optim.residual      =worktn[4,:]
    optim.d             =worktn[5,:]
    optim.Hd            =worktn[6,:]
    optim.eisenvect     =worktn[7,:]
    optim.descent_rec   =worktn[8,:]
    optim.residual_rec  =worktn[9,:]
    #optim.conv_CG=False if optim.conv_CG==0 else True
    if optim.CG_phase!='INIT':
        optim.Hd[:] =np.fromfile(optim.fhess, np.float32)
    
def save_workinfo(foptim, optim):
    optim.worktn.tofile(os.path.join(optim.fdir, 'worktn.bin'))
    optim.d.tofile(optim.fvctr)
    optim.conv_CG=0 if optim.conv_CG==False else 1
    del optim.worktn
    del optim.mod
    del optim.grad
    del optim.descent
    del optim.descent_prev
    del optim.residual
    del optim.d
    del optim.Hd
    del optim.eisenvect
    del optim.descent_rec
    del optim.residual_rec
    del optim.debug
    #save optim
    #print(optim.__dict__)
    #a=optim.__dict__
    #for k in a:
        #print(type(a[k]), k, a[k])
    fout =open(foptim,'w')
    fout.write(json.dumps(optim.__dict__,indent=4))
    fout.close()

def save_final(foptim, optim):
    fdesc=os.path.join(optim.fdir, 'iter'+str(optim.iterc), 'desc.bin')
    optim.descent.tofile(fdesc)
    save_workinfo(foptim, optim)
    print("cls 0: found the descent direction and exit "+ __file__)
    sys.exit(0)
    

#Input paras
fdir    =sys.argv[1] #'optim'
foptim  =os.path.join(fdir,'optim.json')
paras   =eval(open(foptim).read())
paras['fdir']=fdir
#generate a instance of class for work
optim=Dict2Class(paras)
load_worktn(os.path.join(optim.fdir, 'worktn.bin'), optim)

#run linear conjugate gradient algorithm
optim.init_descent = 0
optim.debug=True
if optim.CG_phase=='INIT': #initialize the conjugate gradient process
    optim.residual[:]=optim.grad
    optim.d[:]=-1.*optim.residual
    optim.Hd[:]=0.
    optim.descent[:]=0. #given initial solution
    optim.qk_CG=0. #debug
    optim.hessian_term=0. #debug
    optim.norm_residual =np.linalg.norm(optim.residual).item()
    optim.conv_CG=False
    optim.cpt_iter_CG=0
    optim.CG_phase='IRUN'
    optim.Qp_old=0.0
    optim.min_residual=1.0
    if optim.flag_min_res == 1:
        optim.descent_rec[:]=optim.d
        optim.residual_rec[:]=optim.residual
    if optim.iterc>optim.iter0:
        eta0=optim.eta
        forcing_term_TRN(optim)
        print("cls 0: update eta: ",eta0, optim.eta)
else:#enter linear CG loop
    dHd = np.inner(optim.d,optim.Hd)
    if optim.flag_singularity == 1:
        test_singularity(dHd,optim) #modified by Tao
        if optim.conv_CG == True:
            print('cls 0: test_singularity\n')
            save_final(foptim, optim)
    if dHd<0.:
        #-----------------------------------------------------#
        # if dHd < 0, detection of a negative eigenvalue of   #
        # the Hessian operator, stop the process              #        
        #-----------------------------------------------------#     
        if optim.cpt_iter_CG == 0:
            #-----------------------------------------------------#
            # if this is the first iteration, then return the      #
            # opposite of the gradient as descent direction       #
            # (steepest descent direction)                        #
            #-----------------------------------------------------#     
            optim.descent[:]=optim.d
            #-----------------------------------------------------#
            # if the debug option is activated, compute the       #
            # quadratic function minimized during the conjugate   #
            # gradient process (check is this value decresae      #
            # throughout the CG iterations )                      #
            #-----------------------------------------------------#     
            if optim.debug:
                #allocate(mgrad(n))
                mgrad=-1.*optim.grad
                optim.norm_residual = np.linalg.norm(optim.residual).item()
                alpha=(optim.norm_residual**2)/dHd         
                optim.qkm1_CG=optim.qk_CG
                grad_term = np.inner(optim.descent,mgrad)
                optim.hessian_term=optim.hessian_term+(alpha**2)*dHd
                optim.qk_CG=-grad_term+0.5*optim.hessian_term  
                #deallocate(mgrad)
        optim.conv_CG=True
        print('cls 0: Negative curvature\n')
        save_final(foptim, optim)
    else:
        if optim.flag_negative == 1:        #use test 1A'
            test_negative_curv(dHd,optim) #modified by Tao
            if optim.conv_CG == True:
                print('cls 0: test_negative_curv\n')
                save_final(optim)
        #-----------------------------------------------------#
        # if dHd > 0, then perform one conjugate gradient     #
        # iteration                                           #
        #-----------------------------------------------------#     
        #Update descent direction
        optim.norm_residual = np.linalg.norm(optim.residual).item()
        alpha=(optim.norm_residual**2)/dHd
        #Recording
        optim.descent_prev[:]=optim.descent
        optim.descent[:]=optim.descent+alpha*optim.d
        if optim.flag_negative == 2:    #use test 2A
            test_negative_curv2(optim.grad,optim) #modified by Tao
            if optim.conv_CG == True:
                print('cls 0: test_negative_curv2\n')
                save_final(foptim, optim)
        optim.residual[:]=optim.residual+alpha*optim.Hd  
        #Update CG direction
        norm_residual_prev=optim.norm_residual
        optim.norm_residual=np.linalg.norm(optim.residual).item()
        beta=(optim.norm_residual**2)/(norm_residual_prev**2)
        optim.d[:]=-1.*optim.residual+beta*optim.d
        print('cls 0: alpha,beta %f %f\n'%(alpha,beta))
        
        #Recording
        if optim.flag_min_res == 1:
            record_min_residual(optim)
        #Update iteration counter 
        optim.cpt_iter_CG=optim.cpt_iter_CG+1
        #-----------------------------------------------------#
        # if the debug option is activated, compute the       #
        # quadratic function minimized during the conjugate   #
        # gradient process (check is this value decresae      #
        # throughout the CG iterations )                      #
        #-----------------------------------------------------#        
        if optim.debug:
            #allocate(mgrad(n))
            mgrad=-1.*optim.grad              
            optim.qkm1_CG=optim.qk_CG
            grad_term = np.inner(optim.descent,mgrad)
            descent_scal_Hd = np.inner(optim.descent_prev,optim.Hd)
            optim.hessian_term=optim.hessian_term+(alpha**2)*dHd+\
                                2.*alpha*descent_scal_Hd
            optim.qk_CG=-grad_term+0.5*optim.hessian_term
            #deallocate(mgrad)
        
        #-----------------------------------------------------#
        # Check if the Eisenstat stopping critertion is       #
        # satisfied                                           #
        #-----------------------------------------------------#
        if optim.flag_truncation >=1:
            test_truncation(optim.grad,optim)
            if optim.conv_CG == True:
                print('cls 0: test_truncation\n')
                save_final(foptim, optim)
        optim.conv_CG=(optim.norm_residual <= optim.eta*optim.norm_grad \
                or optim.cpt_iter_CG >= optim.niter_max_CG)
        if optim.conv_CG == True:
            print('cls 0: conv_CG=(norm_residual<=eta*norm_grad or cpt_iter_CG>=niter_max_CG) satisfied\n')
            save_final(foptim, optim)
        print("cls 0: ",optim.norm_residual, optim.eta, optim.norm_grad, optim.cpt_iter_CG,optim.niter_max_CG )
        #-----------------------------------------------------#
        # Print information on the current CG iteration       #
        #-----------------------------------------------------#
save_workinfo(foptim, optim)
print("cls 0: Iter_CG     qk       norm_res    norm_res/||gk||")
print("cls 0: %5d %10.4e %10.4e %.4f"%(optim.cpt_iter_CG,optim.qk_CG,\
    optim.norm_residual,optim.norm_residual/optim.norm_grad))
print('cls 0'+'#'*60+'end\n')

print("Finished...", __file__)
