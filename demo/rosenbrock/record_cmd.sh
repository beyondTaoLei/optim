# cd directory
cd /home/tao/Nutstore_Files/works/optim/demo/rosenbrock/prj

# inversion
python ../S0_prepare.py optim_PSTD
python ../S1_optim_PSTD.py optim_PSTD 1 200

python ../S0_prepare.py optim_PNLCG
python ../S1_optim_PNLCG.py optim_PNLCG 1 50

python ../S0_prepare.py optim_LBFGS
python ../S1_optim_LBFGS.py optim_LBFGS 1 50

python ../S0_prepare.py optim_TRN
python ../S1_optim_TRN.py optim_TRN 1 50

# figure
python ../prob/plot_result.py optim_PSTD 200
python ../prob/plot_result.py optim_PNLCG 50
python ../prob/plot_result.py optim_LBFGS 25
python ../prob/plot_result.py optim_TRN 16
python ../prob/plot_misfits_comp.py

# clean all data
#rm -rf FD fig optim_*
