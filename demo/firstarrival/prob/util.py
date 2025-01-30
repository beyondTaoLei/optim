#!/usr/bin/env python3
import os 
import shutil
import numpy as np

def raypath2senstivity(raypath, xvel, yvel, mat2d):
    mat2d *=0.0
    npoints = len(raypath)
    dx = xvel[1] - xvel[0]
    dy = yvel[1] - yvel[0]
    col, ray = [], []
    for i in range(npoints-1):
        #print(i)
        a = [raypath[i, 0], raypath[i, 1]]
        b = [raypath[i+1, 0], raypath[i+1, 1]]
        x0 = (a[0] + b[0])/2.0
        y0 = (a[1] + b[1])/2.0
        ix = int((x0-xvel[0])//dx)
        iy = int((y0-yvel[0])//dy)
        dist = np.sqrt((a[0]-b[0])**2+(a[1]-b[1])**2)
        mat2d[ix,iy] = dist

def next_xy(gx, x0, xvel):
    """
    the grid coordinate in x axis
    """
    if gx > 0:
        try:
            x = xvel[xvel<x0][-1]
        except IndexError:
            tmp = np.array([x0])
            x = np.clip(tmp, xvel[0], xvel[-1])[0]
    elif gx < 0:
        try:
            x = xvel[xvel>x0][0]
        except IndexError:
            tmp = np.array([x0])
            x = np.clip(tmp, xvel[0], xvel[-1])[0]
    else:
        x = x0
    return x
        

def intersec(x0, y0, gx, gy, xvel, yvel, v_current):
    """
    """
    # find next x, then y
    x = next_xy(gx, x0, xvel)
    d_ray = (x0 - x) / (gx * v_current)
    y = y0- gy * v_current * d_ray
    y = np.clip(np.array([y]), yvel[0], yvel[-1])[0]
    distx = np.sqrt((x0-x)**2+(y0-y)**2)
    resx = [x, y]
    
    # find next y, then x
    y = next_xy(gy, y0, yvel)
    d_ray = (y0 - y) / (gy * v_current)
    x = x0- gx * v_current * d_ray
    x = np.clip(np.array([x]), xvel[0], xvel[-1])[0]
    disty = np.sqrt((x0-x)**2+(y0-y)**2)
    resy = [x, y]
    
    if distx < disty:
        if distx > 0.0:
            return resx
        else:
            return resy
    else:
        if disty > 0.0:
            return resy
        else:
            return resx

def linear_interpolation_2d(xvel, yvel, vmat, x0, y0):
    x0 = np.clip(np.array([x0]), xvel[0], xvel[-1])[0]
    y0 = np.clip(np.array([y0]), yvel[0], yvel[-1])[0]
    dx, dy = xvel[1] - xvel[0], yvel[1] - yvel[0]
    ix1, ix2 = np.argsort(np.abs(x0-xvel))[:2]
    iy1, iy2 = np.argsort(np.abs(y0-yvel))[:2]
    v1 = vmat[ix1,iy1] +(vmat[ix2,iy1] - vmat[ix1,iy1]) / dx * (x0-xvel[ix1])
    v2 = vmat[ix1,iy2] +(vmat[ix2,iy2] - vmat[ix1,iy2]) / dx * (x0-xvel[ix1])
    v0 = v1 + (v2 - v1) / dy * (y0 - yvel[iy1])
    return v0
    
def not_near_source(x_next, y_next, x_src, y_src, d_ray):
    offx = x_next - x_src
    offy = y_next - y_src
    dist = np.sqrt(offx**2 + offy**2)
    return dist > d_ray
    
def FSM_ray(xvel, yvel, x_src, y_src, x_rec, y_rec, vel2d, time2d, gradx, grady, d_ray):
    """
    find the ray path for receiver (x_rec, y_rec)
    copied from Xingzhong Li's code
    """
    x_next = x_rec
    y_next = y_rec
    coords = []
    while not_near_source(x_next, y_next, x_src, y_src, d_ray):
        x_current = x_next
        y_current = y_next
        coords.append([x_current, y_current])
        v_current = linear_interpolation_2d(xvel, yvel, vel2d, x_current, y_current)
        ttime_current = linear_interpolation_2d(xvel, yvel, time2d, x_current, y_current)
        gradx_current = linear_interpolation_2d(xvel, yvel, gradx, x_current, y_current)
        grady_current = linear_interpolation_2d(xvel, yvel, grady, x_current, y_current)
        ##############
        #print("%.2f %.2f %.5e %.5e"%(x_current, y_current, gradx_current, grady_current))
        x_next, y_next = intersec(x_current, y_current, gradx_current, grady_current, \
                                  xvel, yvel, v_current)
        #x_next = x_current - gradx_current * v_current * d_ray
        #y_next = y_current - grady_current * v_current * d_ray
        #ttime_next = linear_interpolation_2d(xvel, yvel, time2d, x_next, y_next)
        
        #d = d_ray
        #while (x_next < xvel[0] or x_next > xvel[-1] or \
            #y_next < yvel[0] or y_next > yvel[-1] or ttime_next > ttime_current):
            #d = 0.5 * d
            #if d < 0.01 * d_ray:
                #break
            #x_next = x_current - gradx_current * v_current * d
            #y_next = y_current - grady_current * v_current * d
            #ttime_next = linear_interpolation_2d(xvel, yvel, time2d, x_next, y_next)
        ##############
    coords.append([x_src, y_src])    
    return np.array(coords, np.float32)

def FSM_rays(fout, xvel, yvel, x_src, y_src, coordR, vel2d, time2d, gradx, grady, d_ray):
    ntr, _ = coordR.shape
    dict1={}
    for itr in range(ntr):
        x_rec, y_rec = coordR[itr,0], coordR[itr,1]
        #print(itr, x_src, y_src, x_rec, y_rec)
        coords = FSM_ray(xvel, yvel, x_src, y_src, x_rec, y_rec, vel2d, time2d, gradx, grady, d_ray)
        dict1['%04d'%itr] = coords
    #save to file
    fout = open(fout, 'w')
    for itr in range(ntr):
        coords = dict1['%04d'%itr]
        npoints = coords.shape[0]
        a = fout.write("%d "%(npoints))
    a = fout.write("\n")
    for itr in range(ntr):
        coords = dict1['%04d'%itr]
        npoints = coords.shape[0]
        for i in range(npoints):
            a = fout.write("%12.3f %12.3f\n"%(coords[i,0], coords[i,1]))
    fout.close()

def delete_data(dir1, stations):
    if os.path.exists(dir1):
        shutil.rmtree(dir1, ignore_errors=True)#folder and *.sac
    os.makedirs(dir1,exist_ok=True)
    for i in range(len(stations)):
        os.makedirs(os.path.join(dir1, stations[i]),exist_ok=True)

def solv_loc(told, aa, bb, ss, dx2, dy2):
    tnow=-1.0*np.ones(3)
    #case 1
    a = 1.0/dx2 + 1.0/dy2
    b = -2.0*aa/dx2 -2.0*bb/dy2
    c = aa*aa/dx2 + bb*bb/dy2- ss
    delta = b*b - 4.0*a*c
    if delta >= 0.0:
        sqrtdelta, ta = np.sqrt(delta), 2.0*a
        tmp1 = (-b-sqrtdelta)/ta
        tmp2 = (-b+sqrtdelta)/ta
        if tmp2 >= aa and tmp2 >= bb:
            tnow[0]= tmp2
        if tmp1 >= aa and tmp1 >= bb:
            tnow[0]= tmp1
    #case 2
    a = 1.0/dx2
    b = -2.0*aa/dx2 
    c = aa*aa/dx2 - ss 
    delta = b*b - 4.0*a*c
    if delta >= 0.0:
        sqrtdelta, ta = np.sqrt(delta), 2.0*a
        tmp1 = (-b-sqrtdelta)/ta
        tmp2 = (-b+sqrtdelta)/ta
        if tmp2 >= aa and tmp2 <= bb:
            tnow[1] = tmp2
        if tmp1 >= aa and tmp1 <= bb:
            tnow[1] = tmp1
    #case 3
    a = 1.0/dy2
    b = -2.0*bb/dy2 
    c = bb*bb/dy2 - ss 
    delta = b*b - 4.0*a*c
    if delta >= 0.0:
        sqrtdelta, ta = np.sqrt(delta), 2.0*a
        tmp1 = (-b-sqrtdelta)/ta
        tmp2 = (-b+sqrtdelta)/ta
        if tmp2 <= aa and tmp2 >= bb:
            tnow[2] = tmp2
        if tmp1 <= aa and tmp1 >= bb:
            tnow[2] = tmp1
    aa=tnow[tnow>=0.0].tolist()
    aa.append(told)
    return min(aa)

def fastsweep(ix1, ix2, iy1, iy2, nx, ny, dx, dy, slow, timetable, source):
    dx2, dy2 = dx*dx, dy*dy
    gap = 0.0
    incx=1 if ix1<ix2 else -1
    incy=1 if iy1<iy2 else -1
    for i in range(ix1, ix2+incx, incx):
      for j in range(iy1, iy2+incy, incy):
        if source[i,j]<0.5:
            ss = slow[i,j]*slow[i,j]
            if i>0 and i<nx-1:
                aa  = min(timetable[i-1,j],timetable[i+1,j])
            elif i==0:
                aa = timetable[i+1,j]
            else:
                aa = timetable[i-1,j]
            if j>0 and j<ny-1:
                bb  = min(timetable[i,j-1],timetable[i,j+1])
            elif j==0:
                bb = timetable[i,j+1]
            else:
                bb = timetable[i,j-1]
            tmp=solv_loc(timetable[i,j], aa, bb, ss, dx2, dy2)
            gap = gap + (timetable[i,j]-tmp)**2
            timetable[i,j] = tmp
        else:
            pass
            #print('here is the source ', i,j,timetable[i,j])
    gap = np.sqrt(gap*dx*dy)
    return gap

def find_nearest(arr, value):
    #https://stackoverflow.com/questions/3843017/efficiently-detect-sign-changes-in-python
    idxs=np.where(np.diff(np.sign(value-arr)))[0]
    if len(idxs)>0:
        label, i0, dist=1, idxs[-1], value-arr[idxs[-1]]
    else:
        label, i0, dist=0, 0, 1.0
    return label, i0, dist

def sourceInitiation(sx, sy, xvel, yvel,timetable, source, slow, ixy):
    nx, ny=timetable.shape
    labelx, i1, distx=find_nearest(xvel, sx)
    labely, j1, disty=find_nearest(yvel, sy)
    if labelx*labely==0:
        raise ValueError('I cannot locate the source in the domain: ', sx, sy)
    i1=i1+1 if distx/(xvel[1]-xvel[0])>0.5 else i1
    j1=j1+1 if disty/(yvel[1]-yvel[0])>0.5 else j1
    ii1, ii2 = max(0,i1-ixy), min(nx,i1+ixy+1)
    jj1, jj2 = max(0,j1-ixy), min(ny,j1+ixy+1)
    for i in range(ii1, ii2):
        for j in range(jj1, jj2):
            timetable[i, j]=np.sqrt((sx-xvel[i])**2+(sy-yvel[j])**2)*slow[i1, j1]
            source[i, j]=1.0

def extract_record(coordR, xvel, yvel,timetable):
    ntr, _ = coordR.shape
    time=np.zeros(ntr, np.float64)
    dx, dy=xvel[1]-xvel[0], yvel[1]-yvel[0]
    for itr in range(ntr):
        labelx, i1, distx=find_nearest(xvel, coordR[itr,0])
        labely, j1, disty=find_nearest(yvel, coordR[itr,1])
        if labelx*labely==0:
            time[itr]=-100.0
            print('I cannot locate the receiver in the domain: ', coordR[itr,0], coordR[itr,1])
        else:
            #print('receiver:', coordR[itr,0], coordR[itr,1], i1, j1)
            wx, wy=distx/dx, disty/dy
            wt1, wt2 = (1.0-wx)*(1.0-wy), (1.0-wx)*wy
            wt3, wt4 = wx*wy, wx*(1.0-wy)
            #print(itr, wt1, wt2, wt3, wt4)
            i2, j2 = i1 + 1, j1 + 1
            time[itr] = wt1*timetable[i1,j1]+wt2*timetable[i1,j2]+ \
                wt3*timetable[i2,j2]+wt4*timetable[i2,j1]
    return time

def linefitting(Xs, Ys):
  '''
  Least square method line fitting Y = a*X + b
  '''
  N = len(Xs)
  l1 = np.sum(Xs)
  l2 = np.sum(Ys)
  l3 = np.sum(Xs*Ys)
  l4 = np.sum(Xs*Xs)
  a = (N * l3 - l1 * l2) / (N * l4 - l1 * l1)
  b = (l4* l2 - l1 * l3) / (N * l4 - l1 * l1)
  ave = l2 / N
  return a, b, ave

def fun(a,b):
    return a+b