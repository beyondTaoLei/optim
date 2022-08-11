#common
workdir=/data/ess-leit/fwi_marmousi2
mydemo=/work/ess-leit/optim/demo/fwi
mycase=/work/ess-leit/optim/demo/fwi/case1
mpirunpy=/work/ess-leit/seis_software/anaconda3/envs/nn/bin/mpirun

# #at PC
# cd /home/tao/Nutstore_Files/works/optim/demo/fwi
# scp -r case1 prob  ess-leit@172.18.6.175:${mydemo}/
# 
# #at HPC server
# # make observed data
# cd ${workdir}
# cp model/true/true.vp model/inv/inv.vp
# cp model/true/true.rho model/inv/inv.rho
# python ${mycase}/S0_prepare_files.py
# python ${mycase}/S0_prepare.py optim
# python ${mydemo}/prob/wavelet.py list/eventinfo.csv 1 6500 0.001 source
# ${mpirunpy} -np 1 python ${mydemo}/prob/filtering.py optim/ list/proc_wav.csv 0
# python ${mycase}/S0_generate_data.py optim

# run fwi
#############
## PSTD test
cp model/init/init.vp model/inv/inv.vp
cp model/init/init.rho model/inv/inv.rho
#1st stage
python ${mycase}/S0_prepare.py pstd
fc=2.5
${mpirunpy} -np 10 python ${mydemo}/prob/filtering.py pstd list/proc_wav.csv ${fc}
${mpirunpy} -np 10 python ${mydemo}/prob/filtering.py pstd list/proc_data.csv ${fc}
python ${mycase}/S1_optim_PSTD.py pstd 1 30

#2nd stage
python ${mycase}/S0_prepare_fc5.py pstd
fc=5.0
${mpirunpy} -np 10 python ${mydemo}/prob/filtering.py pstd list/proc_wav.csv ${fc}
${mpirunpy} -np 10 python ${mydemo}/prob/filtering.py pstd list/proc_data.csv ${fc}
python ${mycase}/S1_optim_PSTD.py pstd 31 60

#3rd stage
python ${mycase}/S0_prepare_fc7.py pstd
fc=7.0
${mpirunpy} -np 10 python ${mydemo}/prob/filtering.py pstd list/proc_wav.csv ${fc}
${mpirunpy} -np 10 python ${mydemo}/prob/filtering.py pstd list/proc_data.csv ${fc}
python ${mycase}/S1_optim_PSTD.py pstd 61 90

#############
## PNLCG test
cp model/init/init.vp model/inv/inv.vp
cp model/init/init.rho model/inv/inv.rho
#1st stage
python ${mycase}/S0_prepare.py pnlcg
fc=2.5
${mpirunpy} -np 10 python ${mydemo}/prob/filtering.py pnlcg list/proc_wav.csv ${fc}
${mpirunpy} -np 10 python ${mydemo}/prob/filtering.py pnlcg list/proc_data.csv ${fc}
python ${mycase}/S1_optim_PNLCG.py pnlcg 1 30

#2nd stage
python ${mycase}/S0_prepare_fc5.py pnlcg
fc=5.0
${mpirunpy} -np 10 python ${mydemo}/prob/filtering.py pnlcg list/proc_wav.csv ${fc}
${mpirunpy} -np 10 python ${mydemo}/prob/filtering.py pnlcg list/proc_data.csv ${fc}
python ${mycase}/S1_optim_PNLCG.py pnlcg 31 60

#3rd stage
python ${mycase}/S0_prepare_fc7.py pnlcg
fc=7.0
${mpirunpy} -np 10 python ${mydemo}/prob/filtering.py pnlcg list/proc_wav.csv ${fc}
${mpirunpy} -np 10 python ${mydemo}/prob/filtering.py pnlcg list/proc_data.csv ${fc}
python ${mycase}/S1_optim_PNLCG.py pnlcg 61 90

#############
## l-BFGS test
cp model/init/init.vp model/inv/inv.vp
cp model/init/init.rho model/inv/inv.rho
#1st stage
python ${mycase}/S0_prepare.py lbfgs
fc=2.5
${mpirunpy} -np 10 python ${mydemo}/prob/filtering.py lbfgs list/proc_wav.csv ${fc}
${mpirunpy} -np 10 python ${mydemo}/prob/filtering.py lbfgs list/proc_data.csv ${fc}
python ${mycase}/S1_optim_LBFGS.py lbfgs 1 30

#2nd stage
python ${mycase}/S0_prepare_fc5.py lbfgs
fc=5.0
${mpirunpy} -np 10 python ${mydemo}/prob/filtering.py lbfgs list/proc_wav.csv ${fc}
${mpirunpy} -np 10 python ${mydemo}/prob/filtering.py lbfgs list/proc_data.csv ${fc}
python ${mycase}/S1_optim_LBFGS.py lbfgs 31 60

#3rd stage
python ${mycase}/S0_prepare_fc7.py lbfgs
fc=7.0
${mpirunpy} -np 10 python ${mydemo}/prob/filtering.py lbfgs list/proc_wav.csv ${fc}
${mpirunpy} -np 10 python ${mydemo}/prob/filtering.py lbfgs list/proc_data.csv ${fc}
python ${mycase}/S1_optim_LBFGS.py lbfgs 61 90

#############
## TRN test
cp model/init/init.vp model/inv/inv.vp
cp model/init/init.rho model/inv/inv.rho
#1st stage
python ${mycase}/S0_prepare.py trn
fc=2.5
${mpirunpy} -np 10 python ${mydemo}/prob/filtering.py trn list/proc_wav.csv ${fc}
${mpirunpy} -np 10 python ${mydemo}/prob/filtering.py trn list/proc_data.csv ${fc}
python ${mycase}/S1_optim_TRN.py trn 1 30

#2nd stage
python ${mycase}/S0_prepare_fc5.py trn
fc=5.0
${mpirunpy} -np 10 python ${mydemo}/prob/filtering.py trn list/proc_wav.csv ${fc}
${mpirunpy} -np 10 python ${mydemo}/prob/filtering.py trn list/proc_data.csv ${fc}
python ${mycase}/S1_optim_TRN.py trn 31 60

#3rd stage
python ${mycase}/S0_prepare_fc7.py trn
fc=7.0
${mpirunpy} -np 10 python ${mydemo}/prob/filtering.py trn list/proc_wav.csv ${fc}
${mpirunpy} -np 10 python ${mydemo}/prob/filtering.py trn list/proc_data.csv ${fc}
python ${mycase}/S1_optim_TRN.py trn 61 90







