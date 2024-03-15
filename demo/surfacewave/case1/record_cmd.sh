conda activate disba
# make observed data
mycase=/home/tao/seis_software/test/optim/demo/surfacewave/case1
cd /home/tao/seis_software/test/optim/demo/surfacewave/prj

cp model/true/true.vp model/inv/inv.vp
cp model/true/true.vs model/inv/inv.vs
cp model/true/true.rho model/inv/inv.rho
python ${mycase}/T0_prepare_files.py
python ${mycase}/S0_prepare.py optim
python ${mycase}/T0_generate_data.py optim

# run tomography
## LSQR test
cp model/init/init.vp model/inv/inv.vp
cp model/init/init.vs model/inv/inv.vs
cp model/init/init.rho model/inv/inv.rho
python ${mycase}/S0_prepare.py lsqr
python ${mycase}/S1_optim_LSQR.py lsqr 1 20
python ${mycase}/S2_run.py lsqr
# python ../prob/plot_L_curve.py lsqr/iter1/L_curve.dat
# ximage n1=400 <lsqr/iter1/L_curve_12.bin cmap=hsv2&

## PLCG test
cp model/init/init.vp model/inv/inv.vp
cp model/init/init.vs model/inv/inv.vs
cp model/init/init.rho model/inv/inv.rho
python ${mycase}/S0_prepare.py plcg
python ${mycase}/S1_optim_PLCG.py plcg 1 20
python ${mycase}/S2_run.py plcg

## STD test
cp model/init/init.vp model/inv/inv.vp
cp model/init/init.vs model/inv/inv.vs
cp model/init/init.rho model/inv/inv.rho
python ${mycase}/S0_prepareNL.py std
python ${mycase}/S1_optim_STD.py std 1 50
python ${mycase}/S2_run.py std

## NLCG test
cp model/init/init.vp model/inv/inv.vp
cp model/init/init.vs model/inv/inv.vs
cp model/init/init.rho model/inv/inv.rho
python ${mycase}/S0_prepareNL.py nlcg
python ${mycase}/S1_optim_NLCG.py nlcg 1 50
python ${mycase}/S2_run.py nlcg

## l-BFGS test
cp model/init/init.vp model/inv/inv.vp
cp model/init/init.vs model/inv/inv.vs
cp model/init/init.rho model/inv/inv.rho
python ${mycase}/S0_prepareNL.py lbfgs
python ${mycase}/S1_optim_LBFGS.py lbfgs 1 30
python ${mycase}/S2_run.py lbfgs

## TRN test
cp model/init/init.vp model/inv/inv.vp
cp model/init/init.vs model/inv/inv.vs
cp model/init/init.rho model/inv/inv.rho
python ${mycase}/S0_prepareNL.py trn
python ${mycase}/S1_optim_TRN.py trn 1 20
python ${mycase}/S2_run.py trn

#phase velocity at iteration 20
cp model/init/init.vp model/inv/inv.vp
cp model/init/init.vs model/inv/inv.vs
cp model/init/init.rho model/inv/inv.rho
mpirun -np 8 python3 ${mycase}/../prob/forward.py lsqr
cp -r syn/vmap_M0.bin vmap_M0_init.bin


cp lsqr/iter20/modg.bin model/inv/inv.vs
mpirun -np 8 python3 ${mycase}/../prob/forward.py lsqr
cp -r syn/vmap_M0.bin vmap_M0_lsqr20.bin

cp plcg/iter20/modg.bin model/inv/inv.vs
mpirun -np 8 python3 ${mycase}/../prob/forward.py plcg
cp -r syn/vmap_M0.bin vmap_M0_plcg20.bin

cp std/iter20/modg.bin model/inv/inv.vs
mpirun -np 8 python3 ${mycase}/../prob/forward.py std
cp -r syn/vmap_M0.bin vmap_M0_std20.bin

cp nlcg/iter20/modg.bin model/inv/inv.vs
mpirun -np 8 python3 ${mycase}/../prob/forward.py nlcg
cp -r syn/vmap_M0.bin vmap_M0_nlcg20.bin

cp lbfgs/iter20/modg.bin model/inv/inv.vs
mpirun -np 8 python3 ${mycase}/../prob/forward.py lbfgs
cp -r syn/vmap_M0.bin vmap_M0_lbfgs20.bin

cp trn/iter20/modg.bin model/inv/inv.vs
mpirun -np 8 python3 ${mycase}/../prob/forward.py trn
cp -r syn/vmap_M0.bin vmap_M0_trn20.bin




# plot
ximage n1=400 cmap=hsv2 bclip=4200 wclip=3000 d1=0.5 d2=0.5 legend=1 hbox=400 wbox=1000 <model/true/true.vs title=true |gnome-screenshot -w -d 1 -B -f mod_true.png
ximage n1=400 cmap=hsv2 bclip=4200 wclip=3000 d1=0.5 d2=0.5 legend=1 hbox=400 wbox=1000 <model/init/init.vs title=initial |gnome-screenshot -w -d 1 -B -f mod_init.png
# plot for LSQR
ximage n1=400 cmap=hsv2 bclip=4200 wclip=3000 d1=0.5 d2=0.5 legend=1 hbox=400 wbox=1000 <lsqr/iter20/modg.bin title="lsqr iter20" |gnome-screenshot -w -d 1 -B -f mod_lsqr20.png
# plot for plcg
ximage n1=400 cmap=hsv2 bclip=4200 wclip=3000 d1=0.5 d2=0.5 legend=1 hbox=400 wbox=1000 <plcg/iter20/modg.bin title="plcg iter20" |gnome-screenshot -w -d 1 -B -f mod_plcg20.png

# plot for std
ximage n1=400 cmap=hsv2 bclip=4200 wclip=3000 d1=0.5 d2=0.5 legend=1 hbox=400 wbox=1000 <std/iter20/mod.bin title="std iter20" |gnome-screenshot -w -d 1 -B -f mod_std20.png

# plot for nlcg
ximage n1=400 cmap=hsv2 bclip=4200 wclip=3000 d1=0.5 d2=0.5 legend=1 hbox=400 wbox=1000 <nlcg/iter20/mod.bin title="nlcg iter20" |gnome-screenshot -w -d 1 -B -f mod_nlcg20.png

# plot for lbfgs
ximage n1=400 cmap=hsv2 bclip=4200 wclip=3000 d1=0.5 d2=0.5 legend=1 hbox=400 wbox=1000 <lbfgs/iter20/mod.bin title="lbfgs iter20" |gnome-screenshot -w -d 1 -B -f mod_lbfgs20.png

# plot for trn
ximage n1=400 cmap=hsv2 bclip=4200 wclip=3000 d1=0.5 d2=0.5 legend=1 hbox=400 wbox=1000 <trn/iter20/mod.bin title="trn iter20" |gnome-screenshot -w -d 1 -B -f mod_trn20.png



ximage n1=29 cmap=hsv2 d1=2 d2=0.5 f1=3 legend=1 hbox=400 wbox=1000 <syn/vmap_M0.bin title="init phase velocity" bclip=3800 wclip=2700 |gnome-screenshot -w -d 1 -B -f phv_init.png
ximage n1=29 cmap=hsv2 d1=2 d2=0.5 f1=3 legend=1 hbox=400 wbox=1000 <syn/vmap_M0.bin title="inverted phase velocity" bclip=3800 wclip=2700 &
ximage n1=29 cmap=hsv2 d1=2 d2=0.5 f1=3 legend=1 hbox=400 wbox=1000 <phv_diff.bin title="diff phase velocity" bclip=50 wclip=-50 &

