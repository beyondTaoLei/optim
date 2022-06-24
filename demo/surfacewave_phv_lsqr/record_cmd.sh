conda activate disba
cp model/true/true.vp model/inv/inv.vp
cp model/true/true.vs model/inv/inv.vs
cp model/true/true.rho model/inv/inv.rho
python ../S0_prepare.py optim
python ../S1_generate_data.py optim


cp model/init/init.vp model/inv/inv.vp
cp model/init/init.vs model/inv/inv.vs
cp model/init/init.rho model/inv/inv.rho
python ../S0_prepare.py optim
python ../S1_optim_LSQR.py optim 1 20


ximage n1=400 cmap=hsv2 bclip=4200 wclip=3000 d1=0.5 d2=0.5 legend=1 hbox=400 wbox=1000 title=iter6 <optim/iter6/modg.bin
