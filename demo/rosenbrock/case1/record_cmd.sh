# cd directory
cd /home/tao/seis_software/test/optim/demo/rosenbrock/prj

# inversion
python ../S0_prepare.py optim_STD
python ../S1_optim_STD.py optim_STD 1 200

python ../S0_prepare.py optim_NLCG
python ../S1_optim_NLCG.py optim_NLCG 1 50

python ../S0_prepare.py optim_LBFGS
python ../S1_optim_LBFGS.py optim_LBFGS 1 50

python ../S0_prepare.py optim_TRN
python ../S1_optim_TRN.py optim_TRN 1 50

# figure
python ../plot/plot_result.py optim_STD 200
python ../plot/plot_result.py optim_NLCG 50
python ../plot/plot_result.py optim_LBFGS 25
python ../plot/plot_result.py optim_TRN 16
python ../plot/plot_misfits_comp.py

# clean all data
#rm -rf FD fig optim_*
