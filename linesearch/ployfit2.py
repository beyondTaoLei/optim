#!/usr/bin/env python3
import numpy as np

#/*
# *------------------------------------------------------------------------
# *
# *   Solve linear equation systems by Python3
# *
# *------------------------------------------------------------------------
#*/
#in this case, we will call this function
#epst(1)=0.0
#epst(2)=alpha_lb
#epst(3)=alpha_lb*beta (beta>1.0)
#f(2)<f(1)<=f(3)
#
def ployfit2(epst_old,f_old,alpha_lb,alpha_ub):
    
    AA = np.zeros((3,3))
    bb = np.zeros((3))
    epst = np.zeros((3))
    f = np.zeros((3))
    
    epst_min = np.amin(epst_old)
    epst_max = np.amax(epst_old)
    
    i = np.argmin(epst_old)
    k = np.argmax(epst_old)
    for j in range(3):
        if j != i and j != k:
            break
    
    epst[0] = epst_old[i]
    epst[1] = epst_old[j]
    epst[2] = epst_old[k]
    f[0] = f_old[i]
    f[1] = f_old[j]
    f[2] = f_old[k]
    
    if f[0] == f[1] and f[1] == f[2]:
        print('There are some mistakes in misfit estimation, please check it!!!')
        raise ZeroDivisionError
    
    #ployfit2=alpha_lb
    for i in range(3):
        AA[i,2]=epst[i]*epst[i]
        AA[i,1]=epst[i]
        AA[i,0]=1.0
        bb[i]=f[i]
    ###############
    #solvelin_f(AA, bb, x, 3, 1)
    x = np.linalg.solve(AA, bb)
    opt=-x[1]/(2*x[2])
    ###############
    #write(*,*)'opt=',opt
    
    if x[2]>0.0:
        ployfit2=opt                            #case 1
        if opt <= 0.0:                          #case 5
            ployfit2=epst[0]
        elif opt>=(5.0*epst[2]):                #case 6
            ployfit2=epst[2]
            
    else:
        ployfit2=epst[0]                        #case 2 #case 4
        dis1=abs(epst[0]-opt);
        dis2=abs(epst[1]-opt);
        dis3=abs(epst[2]-opt);
        if dis1 <= dis2 and dis1 <= dis3:       #case 3
            ployfit2=epst[2]
    
    #write(*,*)'123'
    if(ployfit2 < alpha_lb):   ployfit2=alpha_lb
    if(ployfit2 > alpha_ub):   ployfit2=alpha_ub
    
    return ployfit2
    
if __name__=='__main__':
    """
    Supposing misfit function is y = (x-2)^2 -5, 
    the test results are (1, -4),(3, -4), (5, 4), that is:
    
    epst_old = np.array([1,3,5])
    f_old = np.array([-4,-4,4])
    alpha_lb=-10
    alpha_ub=10
    slt=ployfit2(epst_old,f_old,alpha_lb,alpha_ub)
    print(slt)
    
    and the solution is 2.0
    """

    epst_old = np.array([1,3,5])
    f_old = np.array([-4,-4,4])
    alpha_lb=-10
    alpha_ub=10
    slt=ployfit2(epst_old,f_old,alpha_lb,alpha_ub)
    print(slt)
    