#!/usr/bin/env python3
import os
import sys
import glob
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
fiterc = sys.argv[1]

# merge L_curve*.dat to L_curve.csv
L_list = glob.glob(fiterc+'/L_curve*.dat')
df = pd.concat(map(pd.read_csv, L_list),ignore_index=True)
df = df.sort_values(by='lambda', ascending=False)
num = len(df)
df['idx'] = np.arange(num)
fnm = os.path.join(fiterc, 'L_curve.csv')
df.to_csv(fnm, index=False)
# print(df)

# plot L-curve
idx=np.array(df["idx"])
lam=np.array(df["lambda"])
phid=np.array(df["phid"])
phim=np.array(df["phim"])
bnorm=np.array(df["bnorm"])
print("idx     lam         phid        phim        bnorm    lam*phim  lam*phim/phid  phid/bnorm")
for i in range(len(lam)):
    print("%3d %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e"% \
        (i, lam[i], phid[i], phim[i], bnorm[i], lam[i]*phim[i], lam[i]*phim[i]/phid[i], phid[i]/bnorm[i]))
fout = os.path.join(fiterc, 'L_curve.png')
plot_L_curve(phim, phid, lam, fout)

print("Finished...", __file__)