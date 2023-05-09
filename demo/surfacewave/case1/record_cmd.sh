conda activate disba
# make observed data
mydemo=/home/tao/Nutstore_Files/works/optim/demo/surfacewave_phv
mycase=/home/tao/Nutstore_Files/works/optim/demo/surfacewave_phv/case1
cd /home/tao/Nutstore_Files/works/optim/demo/surfacewave_phv/prj
cp model/true/true.vp model/inv/inv.vp
cp model/true/true.vs model/inv/inv.vs
cp model/true/true.rho model/inv/inv.rho
python ${mycase}/S0_prepare_files.py
python ${mycase}/S0_prepare.py optim
python ${mycase}/S0_generate_data.py optim

# run tomography
## LSQR test
cp model/init/init.vp model/inv/inv.vp
cp model/init/init.vs model/inv/inv.vs
cp model/init/init.rho model/inv/inv.rho
python ${mycase}/S0_prepare.py lsqr
python ${mycase}/S1_optim_LSQR.py lsqr 1 20
# python ../prob/plot_L_curve.py lsqr/iter1/L_curve.dat
# ximage n1=400 <lsqr/iter1/L_curve_12.bin cmap=hsv2&

## PLCG test
cp model/init/init.vp model/inv/inv.vp
cp model/init/init.vs model/inv/inv.vs
cp model/init/init.rho model/inv/inv.rho
python ${mycase}/S0_prepare.py plcg
python ${mycase}/S1_optim_PLCG.py plcg 1 20

## PSTD test
cp model/init/init.vp model/inv/inv.vp
cp model/init/init.vs model/inv/inv.vs
cp model/init/init.rho model/inv/inv.rho
python ${mycase}/S0_prepareNL.py pstd
python ${mycase}/S1_optim_PSTD.py pstd 1 50

## PNLCG test
cp model/init/init.vp model/inv/inv.vp
cp model/init/init.vs model/inv/inv.vs
cp model/init/init.rho model/inv/inv.rho
python ${mycase}/S0_prepareNL.py pnlcg
python ${mycase}/S1_optim_PNLCG.py pnlcg 1 50

## l-BFGS test
cp model/init/init.vp model/inv/inv.vp
cp model/init/init.vs model/inv/inv.vs
cp model/init/init.rho model/inv/inv.rho
python ${mycase}/S0_prepareNL.py lbfgs
python ${mycase}/S1_optim_LBFGS.py lbfgs 1 30

## TRN test
cp model/init/init.vp model/inv/inv.vp
cp model/init/init.vs model/inv/inv.vs
cp model/init/init.rho model/inv/inv.rho
python ${mycase}/S0_prepareNL.py trn
python ${mycase}/S1_optim_TRN.py trn 1 20

#phase velocity at iteration 20
cp model/init/init.vp model/inv/inv.vp
cp model/init/init.vs model/inv/inv.vs
cp model/init/init.rho model/inv/inv.rho
mpirun -np 8 python3 ${mydemo}/prob/forward.py lsqr
cp -r syn/vmap_M0.bin vmap_M0_init.bin


cp lsqr/iter20/modg.bin model/inv/inv.vs
mpirun -np 8 python3 ${mydemo}/prob/forward.py lsqr
cp -r syn/vmap_M0.bin vmap_M0_lsqr20.bin

cp plcg/iter20/modg.bin model/inv/inv.vs
mpirun -np 8 python3 ${mydemo}/prob/forward.py plcg
cp -r syn/vmap_M0.bin vmap_M0_plcg20.bin

cp pstd/iter20/modg.bin model/inv/inv.vs
mpirun -np 8 python3 ${mydemo}/prob/forward.py pstd
cp -r syn/vmap_M0.bin vmap_M0_pstd20.bin

cp pnlcg/iter20/modg.bin model/inv/inv.vs
mpirun -np 8 python3 ${mydemo}/prob/forward.py pnlcg
cp -r syn/vmap_M0.bin vmap_M0_pnlcg20.bin

cp lbfgs/iter20/modg.bin model/inv/inv.vs
mpirun -np 8 python3 ${mydemo}/prob/forward.py lbfgs
cp -r syn/vmap_M0.bin vmap_M0_lbfgs20.bin

cp trn/iter20/modg.bin model/inv/inv.vs
mpirun -np 8 python3 ${mydemo}/prob/forward.py trn
cp -r syn/vmap_M0.bin vmap_M0_trn20.bin




# plot
ximage n1=400 cmap=hsv2 bclip=4200 wclip=3000 d1=0.5 d2=0.5 legend=1 hbox=400 wbox=1000 <model/true/true.vs title=true |gnome-screenshot -w -d 1 -B -f mod_true.png
ximage n1=400 cmap=hsv2 bclip=4200 wclip=3000 d1=0.5 d2=0.5 legend=1 hbox=400 wbox=1000 <model/init/init.vs title=initial |gnome-screenshot -w -d 1 -B -f mod_init.png
# plot for LSQR
ximage n1=400 cmap=hsv2 bclip=4200 wclip=3000 d1=0.5 d2=0.5 legend=1 hbox=400 wbox=1000 <lsqr/iter20/modg.bin title="lsqr iter20" |gnome-screenshot -w -d 1 -B -f mod_lsqr20.png
# plot for plcg
ximage n1=400 cmap=hsv2 bclip=4200 wclip=3000 d1=0.5 d2=0.5 legend=1 hbox=400 wbox=1000 <plcg/iter20/modg.bin title="plcg iter20" |gnome-screenshot -w -d 1 -B -f mod_plcg20.png

# plot for pstd
ximage n1=400 cmap=hsv2 bclip=4200 wclip=3000 d1=0.5 d2=0.5 legend=1 hbox=400 wbox=1000 <pstd/iter20/mod.bin title="pstd iter20" |gnome-screenshot -w -d 1 -B -f mod_pstd20.png

# plot for pnlcg
ximage n1=400 cmap=hsv2 bclip=4200 wclip=3000 d1=0.5 d2=0.5 legend=1 hbox=400 wbox=1000 <pnlcg/iter20/mod.bin title="pnlcg iter20" |gnome-screenshot -w -d 1 -B -f mod_pnlcg20.png

# plot for lbfgs
ximage n1=400 cmap=hsv2 bclip=4200 wclip=3000 d1=0.5 d2=0.5 legend=1 hbox=400 wbox=1000 <lbfgs/iter20/mod.bin title="lbfgs iter20" |gnome-screenshot -w -d 1 -B -f mod_lbfgs20.png

# plot for trn
ximage n1=400 cmap=hsv2 bclip=4200 wclip=3000 d1=0.5 d2=0.5 legend=1 hbox=400 wbox=1000 <trn/iter20/mod.bin title="trn iter20" |gnome-screenshot -w -d 1 -B -f mod_trn20.png



ximage n1=29 cmap=hsv2 d1=2 d2=0.5 f1=3 legend=1 hbox=400 wbox=1000 <syn/vmap_M0.bin title="init phase velocity" bclip=3800 wclip=2700 |gnome-screenshot -w -d 1 -B -f phv_init.png
ximage n1=29 cmap=hsv2 d1=2 d2=0.5 f1=3 legend=1 hbox=400 wbox=1000 <syn/vmap_M0.bin title="inverted phase velocity" bclip=3800 wclip=2700 &
ximage n1=29 cmap=hsv2 d1=2 d2=0.5 f1=3 legend=1 hbox=400 wbox=1000 <phv_diff.bin title="diff phase velocity" bclip=50 wclip=-50 &

