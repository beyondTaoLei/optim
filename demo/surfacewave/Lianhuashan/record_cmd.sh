conda activate disba
# make observed data
cd /home/tao/seis_software/test/paper_surf2Dph/DispInv_Lianhuashan
python ~/Nutstore_private/works/software/optim/demo/surfacewave/Lianhuashan/S0_prepare_data.py

## l-BFGS test
cp model/init/init.vp model/inv/inv.vp
cp model/init/init.vs model/inv/inv.vs
cp model/init/init.rho model/inv/inv.rho
python ~/Nutstore_private/works/software/optim/demo/surfacewave/Lianhuashan/S0_prepare.py lbfgs
python ~/Nutstore_private/works/software/optim/demo/surfacewave/Lianhuashan/S1_optim_LBFGS.py lbfgs 1 1
cp -r syn syn_init
python ~/Nutstore_private/works/software/optim/demo/surfacewave/Lianhuashan/S1_optim_LBFGS.py lbfgs 2 30

# plot model for lbfgs
ximage n1=100 cmap=hsv2 bclip=4000 wclip=2000 d1=0.1 d2=0.1 legend=1 hbox=400 wbox=1000 x1end=6.0 <model/init/init.vs title=initial&
ximage n1=100 cmap=hsv2 bclip=4000 wclip=2000 d1=0.1 d2=0.1 legend=1 hbox=400 wbox=1000 x1end=6.0 <lbfgs/iter20/mod.bin title='lbfgs iter20'&

#plot phase velocity
ximage n1=32 cmap=hsv2 d1=0.1 d2=0.2 f1=0.2 legend=1 hbox=400 wbox=1000 <obs/vmap_M0.bin title="observed phv" bclip=3500 wclip=1800&
ximage n1=32 cmap=hsv2 d1=0.1 d2=0.2 f1=0.2 legend=1 hbox=400 wbox=1000 <syn_init/vmap_M0.bin title="initial phv" bclip=3500 wclip=1800&
ximage n1=32 cmap=hsv2 d1=0.1 d2=0.2 f1=0.2 legend=1 hbox=400 wbox=1000 <syn/vmap_M0.bin title="inverted phv" bclip=3500 wclip=1800&

python

import numpy as np
import matplotlib.pylab as plt
obs = np.fromfile('obs/vmap_M0.bin',np.float32)
inv = np.fromfile('syn/vmap_M0.bin',np.float32)
diff = obs - inv
diff.tofile('diff.bin')
plt.plot(diff[diff>-10000.0],'.')
plt.show()
exit()


ximage n1=32 cmap=hsv2 d1=0.1 d2=0.2 f1=0.2 legend=1 hbox=400 wbox=1000 <diff.bin title="diff phv" bclip=60 wclip=-60&
