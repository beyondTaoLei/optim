#---------------------------------------------------------#
#   suppose all commands are submitted in path ${myprj}   #
#---------------------------------------------------------#

# define paths
myprj=/home/tao/Nutstore_Files/works/optim/demo/firstarrival/prj
mydemo=/home/tao/Nutstore_Files/works/optim/demo/firstarrival
mycase=/home/tao/Nutstore_Files/works/optim/demo/firstarrival/case1
# make observed data
cd ${myprj}
cp model/true/true.vp model/inv/inv.vp
cp model/true/true.vs model/inv/inv.vs
cp model/true/true.rho model/inv/inv.rho
python ${mycase}/S0_prepare_files.py
python ${mycase}/S0_prepare.py optim
python ${mycase}/S0_generate_data.py optim

# run tomography
cp model/init/init.vp model/inv/inv.vp
cp model/init/init.vs model/inv/inv.vs
cp model/init/init.rho model/inv/inv.rho
# model weighting
python ${mycase}/generate_model_weights.py optim

## LSQR test
python ${mycase}/S0_prepare.py optim
python ${mycase}/S1_optim_LSQR.py optim 1 20
## Conjugate Gradient test
python ${mycase}/S0_prepare.py optim_cg
python ${mycase}/S1_optim_LSQR.py optim_cg 1 20

# plot
ximage n1=155 <model/true/true.vp cmap=hsv2 bclip=4000 wclip=1500 d1=25 d2=25 legend=1 hbox=400 wbox=1000 title=true&
ximage n1=155 <model/init/init.vp cmap=hsv2 bclip=4000 wclip=1500 d1=25 d2=25 legend=1 hbox=400 wbox=1000 title=initial&
# plot for LSQR
ximage n1=155 <optim/iter20/modg.bin cmap=hsv2 bclip=4000 wclip=1500 d1=25 d2=25 legend=1 hbox=400 wbox=1000 title="lsqr iter20"&
# plot for CG
ximage n1=155 <optim_cg/iter20/modg.bin cmap=hsv2 bclip=4000 wclip=1500 d1=25 d2=25 legend=1 hbox=400 wbox=1000 title="cg iter20"&
