#!/usr/bin/env python3
import numpy as np
from print_info import print_info

#*****************************************************#
#*                 OPTIMIZATION TOOLBOX              *#
#*****************************************************#
#------------------------------------------------------
# This file includes:                                 #  
#                   1)singularity test                #
#                   2)Negative curvature test         #
#                   3)Truncation test                 #
# in order to improve inner loop CG for truncated     #
#          Newton methods based on TNPACK.            #  
#-----------------------------------------------------#

def test_singularity(dHd,optim):
    """
    Singularity test based on the algorithm from Xie Dexuan.
    @article{xie1999efficient,
        title={Efficient implementation of the truncated-Newton algorithm for large-scale chemistry applications},
        author={Xie, Dexuan and Schlick, Tamar},
        journal={SIAM Journal on Optimization},
        volume={10},
        number={1},
        pages={132--154},
        year={1999},
        publisher={SIAM}
    }
    """
    if(abs(dHd) <= optim.pro_sing):
        #print('test_singularity fails: abs(dHd)<=pro_sing %E %E'%(dHd,optim.pro_sing))
        str1 = 'test_singularity fails: abs(dHd)<=pro_sing %E %E'%(dHd,optim.pro_sing)
        print_info(optim.fdesclog, str1)
        optim.flag_CG_test=1
        exit_inner_loop(optim)    
    
def test_negative_curv(dHd,optim):
    """
    Negative curvature test based on the algorithm from Xie Dexuan.
    xie: Determine descent search directions
    flag_negative = 1:  use Test 1A' (the modified negative curvature test).
                A default choice
    flag_negative = 2:  use Test 2A (the strong negative curvature test)

    @article{xie1999efficient,
        title={Efficient implementation of the truncated-Newton algorithm for large-scale chemistry applications},
        author={Xie, Dexuan and Schlick, Tamar},
        journal={SIAM Journal on Optimization},
        volume={10},
        number={1},
        pages={132--154},
        year={1999},
        publisher={SIAM}
    }
    """
    norm_d2 = np.inner(optim.d, optim.d)
    if dHd <= norm_d2*optim.pro_curv:
        #print('test_negative_curv fails: dHd<=norm_d2*pro_curv %E %E %E %E'\
                    #%(dHd,norm_d2*optim.pro_curv,norm_d2,optim.pro_curv))
        str1 = 'test_negative_curv fails: dHd<=norm_d2*pro_curv %E %E %E %E'%(dHd,norm_d2*optim.pro_curv,norm_d2,optim.pro_curv)
        print_info(optim.fdesclog, str1)
        optim.flag_CG_test=2
        exit_inner_loop(optim)        
    
def test_negative_curv2(grad,optim):
    """
    Negative curvature test based on the algorithm from Xie Dexuan.
    xie: Determine descent search directions
    flag_negative = 1:  use Test 1A' (the modified negative curvature test).
                A default choice
    flag_negative = 2:  use Test 2A (the strong negative curvature test)
    """
    gp_new = np.inner(grad, optim.descent)
    gp_old = np.inner(grad, optim.descent_prev)
    if gp_new >= gp_old:
        #print('test_negative_curv2 fails: gp_new >= gp_old %E %E'%(gp_new,gp_old))
        str1 = 'test_negative_curv2 fails: gp_new >= gp_old %E %E'%(gp_new,gp_old)
        print_info(optim.fdesclog, str1)
        optim.flag_CG_test=3
        exit_inner_loop(optim)        

def test_truncation(grad,optim):
    """
    xie: Determine descent search directions
    flag_truncation = 1:  based on residuals of linear equations
    flag_truncation = 2:  based on parabolic function. 
                          A default choice
    """
    if(optim.flag_truncation == 1):      #
        resp = np.inner(optim.residual,optim.descent)
        gp_new = np.inner(grad,optim.descent)
        optim.Qp = 0.5*(resp+gp_new)
        TMP=(optim.cpt_iter_CG+1)*(1.0-(optim.Qp_old/optim.Qp))
        if TMP <= optim.pro_trun:
            #print('test_truncation1 fails: cpt_iter_CG, Qp_old, Qp, TMP, pro_trun%E %E %E %E %E'\
            #%(optim.cpt_iter_CG,optim.Qp_old,optim.Qp,TMP,optim.pro_trun))
            str1 = 'test_truncation1 fails: cpt_iter_CG, Qp_old, Qp, TMP, pro_trun %d %E %E %E %E'%(optim.cpt_iter_CG,optim.Qp_old,optim.Qp,TMP,optim.pro_trun)
            print_info(optim.fdesclog, str1)
            optim.flag_CG_test=5
            exit_inner_loop(optim)
        optim.Qp_old = optim.Qp
    elif optim.flag_truncation == 2:  #
        #norm_grad = np.linalg.norm(grad)
        norm_grad = optim.norm_grad
        ETA=min(norm_grad,optim.pro_trun2/(optim.cpt_iter+1))
        if optim.norm_residual <= ETA*norm_grad:
            #print('test_truncation2 fails: norm_residual,ETA*norm_grad %E %E %E %E'\
                    #%(optim.norm_residual,ETA,norm_grad,ETA*norm_grad))
            str1 = 'test_truncation2 fails: norm_residual,ETA*norm_grad %E %E %E %E'%(optim.norm_residual,ETA,norm_grad,ETA*norm_grad)
            print_info(optim.fdesclog, str1)
            optim.flag_CG_test=6
            exit_inner_loop(optim)    
    
    if optim.cpt_iter_CG >= optim.niter_max_CG:
        #print('test_truncation3 fails: cpt_iter_CG>=niter_max_CG')
        str1 = 'test_truncation3 fails: cpt_iter_CG>=niter_max_CG'
        print_info(optim.fdesclog, str1)
        optim.flag_CG_test=7
        exit_inner_loop(optim)
    
    
def record_min_residual(optim):
    """
    record minimum residual when applying the inner CG agrithm 
    """
    if optim.norm_residual/optim.norm_grad < optim.min_residual:
        #print('record_min_residual works')
        str1 = 'record_min_residual works'
        print_info(optim.fdesclog, str1)
        optim.min_residual=optim.norm_residual/optim.norm_grad
        optim.descent_rec[:]=optim.descent
        optim.residual_rec[:]=optim.residual #be useful to calculate eta

def exit_inner_loop(optim):
    """
    when the inner loop is terminated during generating the Newton direction, 
    one should exit it with the appropriate descent direction.
    """
    optim.init_descent = 0
    ###################
    #if optim.init_descent == 0:
    if optim.cpt_iter_CG == 0:
        optim.descent[:]=optim.d #-grad
    else:
        if optim.flag_CG_test==3:
            optim.descent[:]=optim.descent_prev
    
    if optim.flag_min_res == 1:
        #print('***exit_inner_loop: 2norm_residual/norm_grad,min_residual %E %E'\
                    #%(optim.norm_residual/optim.norm_grad,optim.min_residual))
        str1 = '***exit_inner_loop: 2norm_residual/norm_grad,min_residual %E %E'%(optim.norm_residual/optim.norm_grad,optim.min_residual)
        print_info(optim.fdesclog, str1)
        #replace descent direction by one when residual is minimum. 
        if (optim.norm_residual/optim.norm_grad) > optim.min_residual:
            optim.descent[:]=optim.descent_rec
            optim.residual[:]=optim.residual_rec
            #print('***exit_inner_loop: Replace descent direction by one when residual is minimum: %E'%optim.min_residual)
            str1 = '***exit_inner_loop: Replace descent direction by one when residual is minimum: %E'%optim.min_residual
            print_info(optim.fdesclog, str1)
    optim.conv_CG=True
    #print('***exit_inner_loop: Inner loop is truncated')
    str1 = '***exit_inner_loop: Inner loop is truncated'
    print_info(optim.fdesclog, str1)
    