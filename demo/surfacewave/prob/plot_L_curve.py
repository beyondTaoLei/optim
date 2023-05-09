#!/usr/bin/env python3
import os, sys
import numpy as np
import pandas as pd
import matplotlib.pylab as plt

#########################################################################
# define: plot L-curve 
# usage:
#########################################################################
# functions
def plot_L_curve(phim, phid, lam, fout):
    print("plot L curve ...")
    fig, ax = plt.subplots()
    plt.plot(phim,phid,'r-o')
    plt.grid(True)
    plt.xlabel('x-norm', fontsize=15)
    plt.ylabel('misfit', fontsize=15)
    ax.set_xlim(0, 1.1*max(phim))
    ax.set_title('L-curve', fontsize=17)
    for i in range(len(lam)):
        ax.text(phim[i],phid[i],"(%.2e, %.3f)"%(lam[i], lam[i]*phim[i]/phid[i]))
    plt.savefig(fout)
    plt.show()

#Input paras
fcurve      =sys.argv[1]
fout        =fcurve.replace(fcurve.split('/')[-1],'L_curve.eps')
pd_data=pd.read_csv(fcurve)
idx=np.array(pd_data["idx"])
lam=np.array(pd_data["lambda"])
phid=np.array(pd_data["phid"])
phim=np.array(pd_data["phim"])
bnorm=np.array(pd_data["bnorm"])
print("idx     lam         phid        phim        bnorm    lam*phim  lam*phim/phid  phid/bnorm")
for i in range(len(lam)):
    print("%3d %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e"% \
        (i, lam[i], phid[i], phim[i], bnorm[i], lam[i]*phim[i], lam[i]*phim[i]/phid[i], phid[i]/bnorm[i]))
plot_L_curve(phim, phid, lam, fout)

print("Finished...", __file__)