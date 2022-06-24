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


ximage n1=155 cmap=hsv2 bclip=4000 wclip=1500 d1=25 d2=25 legend=1 hbox=400 wbox=1000 title=iter20 < optim/iter20/modg.bin

ximage n1=155 <model/true/true.vp cmap=hsv2 bclip=4000 wclip=1500 d1=25 d2=25 legend=1 hbox=400 wbox=1000 title=true&
ximage n1=155 <model/init/init.vp cmap=hsv2 bclip=4000 wclip=1500 d1=25 d2=25 legend=1 hbox=400 wbox=1000 title=initial&
ximage n1=155 <optim/iter20/modg.bin cmap=hsv2 bclip=4000 wclip=1500 d1=25 d2=25 legend=1 hbox=400 wbox=1000 title=iter20&
