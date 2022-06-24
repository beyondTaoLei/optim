#figure
python3 prob/plot_result.py optim_PSTD 200
python3 prob/plot_result.py optim_PNLCG 50
python3 prob/plot_result.py optim_LBFGS 25
python3 prob/plot_result.py optim_TRN 16
python3 prob/plot_misfits_comp.py

# clean all data
#rm -rf FD fig optim_*