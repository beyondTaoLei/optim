#!/usr/bin/env python3
import numpy as np
"""
@article{metivier_seiscope_2016,
	title = {The {SEISCOPE} optimization toolbox: A large-scale nonlinear optimization library based on reverse communication},
	volume = {81},
	issn = {0016-8033, 1942-2156},
	url = {http://library.seg.org/doi/10.1190/geo2015-0031.1},
	doi = {10.1190/geo2015-0031.1},
	shorttitle = {The {SEISCOPE} optimization toolbox},
	pages = {F1--F15},
	number = {2},
	journaltitle = {{GEOPHYSICS}},
	author = {Métivier, Ludovic and Brossier, Romain},
	urldate = {2019-01-04},
	date = {2016-03},
	langid = {english}
}
"""
def model_init():
    return np.array([0.25, 0.25], np.float32)

def model_true():
    return np.array([1, 1], np.float32)

def func(xy):
    x,y=xy[0], xy[1]
    return (1-x)**2 +100*(y-x**2)**2

def grad(xy):
    x,y=xy[0], xy[1]
    return np.array([2*(x-1) -400*x*(y-x**2), 200*(y -x**2)], np.float32)

#---------------------------------------------#
#  The routine Rosenbrock_Hess returns        #
#  Hessian-vector product H(x)d in output Hd  #
#  for input parameters x and d               #
#  H is the Hessian matrix                    #
#  x=(x1,x2), d=(d1,d2) are two vector of R^2 #
#---------------------------------------------#
def hess(xy,d):
  #fcost=(1-x)**2+100.*(y-x**2)**2
  #grad[0]=2.*(x-1)-400.*x*(y-x**2)
  #grad[1]=200.*(y-x**2)
  #IN
  #real,dimension[1] :: x,d
  #IN/OUT
  #real,dimension[1] :: Hd
  x,y=xy[0], xy[1]
  Hd1=(1200.*x**2-400.*y+2.)*d[0]-400*x*d[1]
  Hd2=-400.*x*d[0]+200.*d[1]
  return np.array([Hd1,Hd2], np.float32)

def hess2(xy,d):
    """
    @article{schlick_tnpacktruncated_1992,
	title = {{TNPACK}—A Truncated Newton Minimization Package for Large-Scale Problems: I. Algorithm and Usage},
	volume = {18},
	issn = {0098-3500},
	url = {https://doi.org/10.1145/128745.150973},
	doi = {10.1145/128745.150973},
	pages = {46--70},
	number = {1},
	journaltitle = {{ACM} Trans. Math. Softw.},
	author = {Schlick, Tamar and Fogelson, Aaron},
	date = {1992-03}
    }
    """
    #x,y=xy[0], xy[1]
    #MCHPR1 =1.E-17
    epsilon=1.0E-4
    MCHPR1 =2.0*epsilon*(1.0+np.linalg.norm(xy))/np.linalg.norm(d)
    gradp=grad(xy)
    gradq=grad(xy+MCHPR1*d)
    Hd1=(gradq[0]-gradp[0])/MCHPR1
    Hd2=(gradq[1]-gradp[1])/MCHPR1
    print(xy, d,MCHPR1, gradp,gradq)
    return np.array([Hd1,Hd2], np.float32)
