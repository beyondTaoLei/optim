#!/usr/bin/env python3
"""! @brief Manual for programs sofi2Dac and sofi2Dac_grad."""


##
# @mainpage Doxygen Example Project
#
# @section description_main Description
# An example Python program demonstrating how to use Doxygen style comments for
# generating source code documentation with Doxygen.
# https://www.woolseyworkshop.com/2020/06/25/documenting-python-programs-with-doxygen/
#
# @section notes_main Notes
# - Add special project notes here that you want to communicate to the user.
#
# Copyright (c) 2021 ZhangWei's Group (zhangwei@sustech.edu.cn or leit@sustech.edu.cn).  All rights reserved.


##
# @file README.py
#
# @brief Example Python program with Doxygen style comments.
#
# @section description_README Description
# Example Python program with Doxygen style comments.
#
# @section libraries_main Libraries/Modules
# - time standard library (https://docs.python.org/3/library/time.html)
#   - Access to sleep function.
# - sensors module (local)
#   - Access to Sensor and TempSensor classes.
#
# @section notes_README Notes
# - Comments are Doxygen compatible.
#
# @section todo_README TODO
# - None.
#
# @section author_README Author(s)
# - Created by Tao Lei on 05/27/2021.
# - Modified by Tao Lei on 06/11/2021.
#
# Copyright (c) 2021 ZhangWei's Group (zhangwei@sustech.edu.cn or leit@sustech.edu.cn).  All rights reserved.


# Imports
from time import sleep


# Functions
def Areadme():
    """!
    Here you will learn how to calculate seismograms, misfit, adjoint 
    source and gradient for one shot or a given acquisition geometry. 
    Because we run FD with SOFI2DAC, which is modified based on IFOS2D 
    in the acoustic case, you can freely download it from git repository
    https://git.scc.kit.edu/GPIAG-Software/IFOS2D. 
    
    In order to obtain the above results for one shot, you will run 
    the program sofi2Dac, sofi2Dac_grad, calc_misfit_res.py, 
    in which the first two is made in ../src directory and the last 
    one is in the current directory.
    
    # keep in mind, ALL python and bash scripts are submitted in your work directory.
    """
    pass

def creat_project():
    """! Describe how to creat a project for the test.
    Suppose you make a project in directory \b Marmousi2, which includes
    subdirectories as,
        
    @param list     the folder includes \b a. source.dat, receiver.dat, 
    eventname.list and stationname.list, \b b. json files FW.json and GRAD.json, 
    and \b c. the scripts run_seismo.sh, run_misfit.sh, run_misfit_grad.sh for each shot.
    @param model    the folder includes subdirectories as \b true, \b init and \b inv.
    @param source   the folder includes all source time functions.
    @param log      the folder records the program log-file.
    """
    pass
    
def file_MODEL(MODEL = "inv.vp and inv.vs"):
    """! Describe the input model files.
    the models is put in directory <em>model/inv </em> for \b sofi2Dac 
    and \b sofi2Dac_grad. Additionally, we put the true models in 
    <em>model/true </em> with name <em>true.vp</em> and <em>true.rho</em> 
    and put the initial models in <em>model/init </em> with name 
    <em>init.vp</em> and <em>init.rho</em>
    
    @param n1   the grid size in depth direction.
    @param n2   the grid size in horizontal direction.
    @param inv.vp   the vp model file in BIN-format.
    @param inv.rho  the density model file in BIN-format.
    """
    pass
    
def file_SOURCE(SOURCE = "wavelet.su"):
    """! Describe the source time function (STF) file.
    based on previous experience, we simply insert the STF by file for
    different kinds of simulation, as in forward and backward case. here the 
    STF is written in \b SU format, and the position are necessary in the 
    header struct. you can generate several shapes of wavelets with function
    \b write_wavelet.
    """
    pass
    
def file_RECEIVER(RECEIVER = "receiver.dat"):
    """! Describe the receiver position file.
    the receiver positions are given in ASCII format, which is as \n
    500.00	20.00 \n
    520.00	20.00 \n
    540.00	20.00 \n
    560.00	20.00 \n
    580.00	20.00 \n
    600.00	20.00 \n
    620.00	20.00 \n
    ...               \n
    the first column is the \b X coordinate in meter and the second 
    is the \b Y coordinate
    """
    pass

def file_SOURCE_all(SOURCE_all='source.dat'):
    """! Describe all source position file.
    the first row is the number of the source, and the following rows define
    each source function, in which each column is as \n
    coord. \b X, \b tmp, coord. \b Y, time delay \b tshift, peak or main freq. \b f0, amplitude \b amp0.
    
    
    80 \n
    700.00        0.00    20.00    0.20   15.00   1.00e+10 \n
    900.00        0.00    20.00    0.20   15.00   1.00e+10 \n
    1100.00        0.00    20.00    0.20   15.00   1.00e+10 \n
    1300.00        0.00    20.00    0.20   15.00   1.00e+10 \n
    1500.00        0.00    20.00    0.20   15.00   1.00e+10 \n
    1700.00        0.00    20.00    0.20   15.00   1.00e+10 \n
    1900.00        0.00    20.00    0.20   15.00   1.00e+10 \n
    2100.00        0.00    20.00    0.20   15.00   1.00e+10 \n
    2300.00        0.00    20.00    0.20   15.00   1.00e+10 \n
    ...  \n
    
    you can generate all SU-format sources with extern execution write_wavelet, as \n
    
    ${spfroot}/bin/write_wavelet nt=6500 delta=1.0e-03 fsrc=list/source.dat fpath=source
    """
    pass

def file_SOURCENAME(SOURCENAME="eventname.list"):
    """! Describe the source name file.
    Because our propagator is only worked for one shot, we need to distinguish
    the output with different name, the source name file is as \n
    1 sz001 \n
    2 sz002 \n
    3 sz003 \n
    4 sz004 \n
    5 sz005 \n
    6 sz006 \n
    7 sz007 \n
    8 sz008 \n
    9 sz009 \n
   10 sz010 \n
   ...  \n
    in which the first column is the source index, and the second is the source name.
    please keep the same order with \b SOURCE_all.
    
    """
    pass

def file_RECEIVERNAME(RECEIVERNAME="stationname.list"):
    """! Describe the receiver name file.
    the receiver name file is valid for tomography case,which is as \n
    1 wh001 \n
    2 wh002 \n
    3 wh003 \n
    4 wh004 \n
    5 wh005 \n
    6 wh006 \n
    7 wh007 \n
    8 wh008 \n
    9 wh009 \n
   10 wh010 \n
   ...  \n
    in which the first column is the receiver index, and the second is the receiver name.
    please keep the same order with \b RECEIVER.
    """
    pass

def file_FORWARD_JSON(FORWARD_JSON = "FW.json"):
    """! Describe the input json file for forward simulation \b sofi2Dac.
    Because we modify the acoustic simulation from IFOS2d software, the 
    most of input parameters is kept unchanged, which is as follows,\n
    {\n
    "Domain Decomposition" : "comment",\n
                            "NPROCX" : "4",\n
                            "NPROCY" : "1",\n
\n
    "FD order" : "comment", \n
                            "FDORDER" : "8", \n
                            "FDORDER_TIME" : "2", \n
                            "MAXRELERROR" : "1", \n
\n
    "2-D Grid" : "comment", \n
                            "NX" : "1700", \n
                            "NY" : "350", \n
                            "DH" : "10.0", \n
\n
    "Model" : "comment", \n
                            "READMOD" : "1", \n
                            "MFILE" : "../../model/inv/inv", \n
                            "WRITE_MODELFILES" : "0", \n
\n
    "Q-approximation" : "comment", \n
                            "L" : "0", \n
                            "FL1" : "5.0", \n 
                            "TAU" : "0.00001", \n
\n
    "Boundary Conditions" : "comment", \n
                            "FREE_SURF" : "1", \n
                            "BOUNDARY" : "0", \n
\n
                            "FW" : "15", \n
                            "ABS_TYPE" : "1", \n
                            "ABS_TYPE values : CPML-Boundary=1; Damping-Boundary=2" : "comment", \n
\n
    "Parameter for CPML (ABS_TYPE=1)" : "comment", \n
                            "NPOWER" : "4.0", \n
                            "K_MAX_CPML" : "1.0", \n
                            "VPPML" : "3000.0", \n
                            "FPML" : "15.0", \n
\n
    "Parameter for ABS_TYPE=2" : "comment", \n
                            "DAMPING" : "8.0", \n
\n
    "Time Stepping" : "comment", \n
                            "TIME" : "6.50", \n
                            "DT" : "1.0e-03", \n
\n
    "Source" : "comment", \n
                            "FLAG_SU_BIN" : "1", \n
                            "SRC_F_SP" : "./wavelet.su", \n
\n
    "Receiver" : "comment", \n
                            "SEISMO" : "2", \n
                            "READREC" : "1", \n
                            "REC_FILE" : "../comm/receiver.dat",\n
                            "REFRECX, REFRECY" : "0.0 , 0.0",\n
                            "XREC1,YREC1" : "270.0 , 270.0",\n
                            "XREC2,YREC2" : "2000.0 , 270.0",\n
                            "NGEOPH" : "1",\n
                            "SEISMO2" : "0",\n
\n
    "Receiver array" : "comment",\n
                            "REC_ARRAY" : "0",\n
                            "REC_ARRAY_DEPTH" : "70.0",\n
                            "REC_ARRAY_DIST" : "40.0", \n
                            "DRX" : "4",\n
\n
    "Seismograms" : "comment",\n
                            "NDT" : "1",\n
                            "SEIS_FORMAT" : "1",\n
                            "SEIS_FILE" : "./seis",\n
\n
    "Snapshots" : "comment",\n
                            "SNAP" : "0",\n
                            "TSNAP1" : "0.0e-3",\n
                            "TSNAP2" : "6.5",\n
                            "TSNAPINC" : "5.0e-03",\n
                            "IDX" : "1",\n
                            "IDY" : "1",\n
                            "SNAP_FORMAT" : "3",\n
                            "SNAP_FILE" : "./snap",\n
                            "FREQ_FILE" : "./freqs_fd.dat",\n
\n                            
    "General inversion parameters" : "comment",\n
                            "FORWARD_ONLY" : "1",\n
                            "TSHIFT" : "0.2",\n
                            "INV_VP" : "1",\n
                            "INV_RHO" : "0",\n
                            "EPRECOND" : "1",\n
                            "EPSILON" : "0.03",\n
\n
    "Monitoring the simulation" : "comment",\n
                            "LOG_FILE" : "test.log",\n
                            "LOG" : "1",\n
                            "OUT_TIMESTEP_INFO" : "10000",\n
\n
    "Checkpoints" : "comment",\n
                            "CHECKPTREAD" : "0",\n
                            "CHECKPTWRITE" : "0",\n
                            "CHECKPT_FILE" : "tmp/checkpoint_sofi2D"\n
    }\n
    The more details you can read the manual of IFOS2D software
    """
    pass

def file_GRADIENT_JSON(GRADIENT_JSON = "GRAD.json"):
    """! Describe the input json file for gradient calculation \b sofi2Dac_grad.
    Actually, the parameters are same with FORWARD_JSON but \b SEISMO2 
    should switch into \b 1, which tells \b sofi2Dac to keep the source 
    wavefields on the boundary for future gradient calculation.
    """
    pass

def run_FORWARD_SCRIPT_one_shot(FORWARD_SCRIPT = "run_seismo.sh"):
    """! Describe how to run forwad modelling to get seismograms.
    This script \b run_seismo.sh is as \n
    ########################################################### \n
    #!/bin/bash \n
    mpiexe='/share/intel/2018u4/compilers_and_libraries_2018.5.274/linux/mpi/intel64/bin/mpiexec' \n
    spfroot='/work/ess-leit/seisplatform' \n
    NP=4 \n
    rm -f seis_p_shot1.su \n
    #forward modelling\n
    ${mpiexe} -np ${NP} ${spfroot}/bin/sofi2Dac ../comm/FW.json \n
    ########################################################### \n
    @param MODEL The files for \b sofi2Dac: <em>inv.vp, inv.rho</em>, if \b necessary, you can update them.
    @param SOURCE wavelet.su
    @return  seismograms (seis_p_shot1.su)
    """
    pass

def run_MISFIT_SCRIPT_one_shot(MISFIT_SCRIPT = "run_misfit.sh"):
    """! Describe how to calculate the misfit between synthetic and observed data. 
    This script \b run_misfit.sh is as \n
    ########################################################### \n
    #!/bin/bash\n
    mpiexe='/share/intel/2018u4/compilers_and_libraries_2018.5.274/linux/mpi/intel64/bin/mpiexec'\n
    spfroot='/work/ess-leit/seisplatform'\n
    NP=4\n

    rm -f seis_p_shot1.su misfit.json\n
    #forward modelling\n
    ${mpiexe} -np ${NP} ${spfroot}/bin/sofi2Dac ../comm/FW.json\n
    mv seis_p_shot1.su syn.su\n
    #data manipulation in future\n
    cp syn.su syn_final.su\n
    #misfit and no adjoint src\n
    python3 ${spfroot}/FD/SOFI2DAC/pyfiles/calc_misfit_res.py 0\n
    ########################################################### \n
    @param MODEL The files for \b sofi2Dac: <em>inv.vp, inv.rho</em>, if \b necessary, you can update them.
    @param SOURCE wavelet.su
    @param DATA observed data \b obs_p_final.su
    @return  misfit (misfit.json)
    """
    pass

def run_GRADIENT_SCRIPT_one_shot(GRADIENT_SCRIPT = "run_misfit_grad.sh"):
    """! Describe how to calculate the gradient for one shot. 
    This script \b run_misfit_grad.sh integrate the whole commands, 
    as forward modelling, misfit and adjoint source calculation, and 
    backward propagation \n
    ########################################################### \n
    #!/bin/bash \n
    mpiexe='/share/intel/2018u4/compilers_and_libraries_2018.5.274/linux/mpi/intel64/bin/mpiexec'\n
    spfroot='/work/ess-leit/seisplatform'\n
    NP=4\n

    rm -f seis_p_shot1.su grad*.bin precond*.bin misfit.json\n
    #forward modelling \n
    ${mpiexe} -np ${NP} ${spfroot}/bin/sofi2Dac ../comm/GRAD.json\n
    mv seis_p_shot1.su syn.su\n
    #data manipulation in future\n
    cp syn.su syn_final.su\n
    #misfit and adjoint src\n
    python3 ${spfroot}/FD/SOFI2DAC/pyfiles/calc_misfit_res.py 1\n
    #gradient\n
    ${mpiexe} -np ${NP} ${spfroot}/bin/sofi2Dac_grad ../comm/GRAD.json\n
    ########################################################### \n
    @param MODEL The files for \b sofi2Dac and \b sofi2Dac_grad: <em>inv.vp, inv.rho</em>, if \b necessary, you can update them.
    @param SOURCE wavelet.su
    @param DATA observed data obs_p_final.su
    @return  misfit and gradient (misfit.json, gradvp_shot1.bin, and gradrho_shot1.bin)
    """
    pass

def Sa_generate_commands_one_iteration():
    """! Describe how to calculate seismograms for all shots.
    
    Firstly, the following files  should be added.\n
    
    in the \b list directory\n
    
    @param SOURCE_all The source information file: <em>source.dat</em>.
    @param RECEIVER The receiver position file: <em>receiver.dat</em>.
    @param FORWARD_JSON The json input files for \b sofi2Dac program: <em>FW.json</em>.
    @param GRADIENT_JSON The json input files for \b sofi2Dac_grad program: <em>GRAD.json</em>.
    @param FORWARD_SCRIPT The script for seismograms: <em>run_seismo.sh</em>
    @param GRADIENT_SCRIPT The script for misfit and gradient: <em>run_misfit_grad.sh</em>
    @param MISFIT_SCRIPT The script for misfit only: <em>run_misfit.sh</em>.
    
    in the \b model/inv directory\n
    @param MODEL The files for \b sofi2Dac and \b sofi2Dac_grad: <em>inv.vp, inv.vs</em>.
    
    After obtaining all the files above, just go to the script S0_prepare.py, and modify 
    the parameters in python dictory [paras] based on your case. The most
    elements is copied from <em>GRAD.json</em>, and run it in your work directory as 
    
    # python /work/ess-leit/seisplatform/FD/SOFI2DAC/pyfiles/S0_prepare.py
    
    This command will generate all scripts for your project, including
    
    Forward modelling: \b cmds_true.sh
    
    Misfit and gradient: \b cmds_misf_grad.sh
    
    Step length estimation: \b cmds_misf_test0.sh and \b cmds_misf_test.sh
    
    Now I will introduce the above scripts one by one.
    
    \b a) \b cmds_true.sh is used to generate the observed data for all shots, and the 
    main task is calling the script \b run_seismo.sh, the content is as
    
    #group# \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz001 && sh run_seismo.sh && mv seis_p_shot1.su obs.su \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz002 && sh run_seismo.sh && mv seis_p_shot1.su obs.su \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz003 && sh run_seismo.sh && mv seis_p_shot1.su obs.su \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz004 && sh run_seismo.sh && mv seis_p_shot1.su obs.su \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz005 && sh run_seismo.sh && mv seis_p_shot1.su obs.su \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz006 && sh run_seismo.sh && mv seis_p_shot1.su obs.su \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz007 && sh run_seismo.sh && mv seis_p_shot1.su obs.su \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz008 && sh run_seismo.sh && mv seis_p_shot1.su obs.su \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz009 && sh run_seismo.sh && mv seis_p_shot1.su obs.su \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz010 && sh run_seismo.sh && mv seis_p_shot1.su obs.su \n
    #group# \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz011 && sh run_seismo.sh && mv seis_p_shot1.su obs.su \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz012 && sh run_seismo.sh && mv seis_p_shot1.su obs.su \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz013 && sh run_seismo.sh && mv seis_p_shot1.su obs.su \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz014 && sh run_seismo.sh && mv seis_p_shot1.su obs.su \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz015 && sh run_seismo.sh && mv seis_p_shot1.su obs.su \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz016 && sh run_seismo.sh && mv seis_p_shot1.su obs.su \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz017 && sh run_seismo.sh && mv seis_p_shot1.su obs.su \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz018 && sh run_seismo.sh && mv seis_p_shot1.su obs.su \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz019 && sh run_seismo.sh && mv seis_p_shot1.su obs.su \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz020 && sh run_seismo.sh && mv seis_p_shot1.su obs.su \n
    ... \n
    
    \b b) \b cmds_misf_grad.sh is used to generate the misfit and gradient for all shots, and the 
    main task is calling the script \b run_misfit_grad.sh, the content is as
    
    #group# \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz001 && sh run_misfit_grad.sh \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz002 && sh run_misfit_grad.sh \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz003 && sh run_misfit_grad.sh \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz004 && sh run_misfit_grad.sh \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz005 && sh run_misfit_grad.sh \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz006 && sh run_misfit_grad.sh \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz007 && sh run_misfit_grad.sh \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz008 && sh run_misfit_grad.sh \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz009 && sh run_misfit_grad.sh \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz010 && sh run_misfit_grad.sh \n
    #group# \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz011 && sh run_misfit_grad.sh \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz012 && sh run_misfit_grad.sh \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz013 && sh run_misfit_grad.sh \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz014 && sh run_misfit_grad.sh \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz015 && sh run_misfit_grad.sh \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz016 && sh run_misfit_grad.sh \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz017 && sh run_misfit_grad.sh \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz018 && sh run_misfit_grad.sh \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz019 && sh run_misfit_grad.sh \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz020 && sh run_misfit_grad.sh \n
    ... \n
    
    #group# \n
    python /work/ess-leit/seisplatform/FD/SOFI2DAC/pyfiles/Sb4_sum_up_misfit.py all \n
    python /work/ess-leit/seisplatform/FD/SOFI2DAC/pyfiles/Sb5_sum_up_gradient.py \n
    python /work/ess-leit/seisplatform/FD/SOFI2DAC/pyfiles/Sb6_process_gradient.py \n
    python /work/ess-leit/seisplatform/FD/SOFI2DAC/pyfiles/Sb7_scale_gradient.py \n
    ... \n
    
    \b c) \b cmds_misf_test0.sh is used to sum up the misfit from the test shots at current model, 
    and the main task is calling the script \b Sb4_sum_up_misfit.py, the content is as
    
    #group# \n
    python /work/ess-leit/seisplatform/FD/SOFI2DAC/pyfiles/Sb4_sum_up_misfit.py test \n
    
    the test source index is defined by \b TESTSHOTINC in S0_prepare.py
    
    \b d) \b  cmds_misf_test.sh is used to sum up the misfit from the test shots at test model, 
    and the content is as
        
    #group# \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz001 && sh run_misfit.sh \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz004 && sh run_misfit.sh \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz007 && sh run_misfit.sh \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz010 && sh run_misfit.sh \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz013 && sh run_misfit.sh \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz016 && sh run_misfit.sh \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz019 && sh run_misfit.sh \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz022 && sh run_misfit.sh \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz025 && sh run_misfit.sh \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz028 && sh run_misfit.sh \n
    #group# \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz031 && sh run_misfit.sh \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz034 && sh run_misfit.sh \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz037 && sh run_misfit.sh \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz040 && sh run_misfit.sh \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz043 && sh run_misfit.sh \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz046 && sh run_misfit.sh \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz049 && sh run_misfit.sh \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz052 && sh run_misfit.sh \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz055 && sh run_misfit.sh \n
    cd /data/ess-leit/Marmousi2_new/fdio/sz058 && sh run_misfit.sh \n
    ...  \n
    #group# \n
    python /work/ess-leit/seisplatform/FD/SOFI2DAC/pyfiles/Sb4_sum_up_misfit.py test
    
    the test model is offered by \b mod_OPTIM modules.
    
    """
    pass
    
def Sb_split_commands_one_iteration(cmd1="cmds_true.sh", cmd2="cmds_misf_grad.sh", cmd3="cmds_misf_test.sh"):
    """! Describe how to split calculation tasks for one iteration.
    Because both FD modelling and gradient calculation need heavy HPC resources, we need adjust
    sub-jobs into different HPC nodes according to the cluster requirements. the commands are as
    
    # python /work/ess-leit/seisplatform/FD/SOFI2DAC/pyfiles/split_cmd.py cmds_true.sh
    
    # python /work/ess-leit/seisplatform/FD/SOFI2DAC/pyfiles/split_cmd.py cmds_misf_grad.sh
    
    # python /work/ess-leit/seisplatform/FD/SOFI2DAC/pyfiles/split_cmd.py cmds_misf_test.sh
    
    here we split jobs based on the flag \b #group#
    @return the scripts (cmds_true_p*.lsf, cmds_misf_grad_p*.lsf, and cmds_misf_test_p*.lsf)
    
    here * denotes the index of sub-jobs.
    """
    pass

def Sc_run_observed_data(cmd="cmds_true.sh or cmds_true_p*.lsf"):
    """! Describe how to generate the observed data for FWI.
    In order to generate the observed data, you should prepare the info.
    as \b model, the source time function \b wavelet.su, the receiver 
    position file, and the input json file for forward simulation \b FW.json.
    the whole steps are:
    
    \b a) copy your true models into model/inv/, which is as
    
    cd model/ \n
    cp true/true.vp inv/inv.vp \n
    cp true/true.vs inv/inv.vs \n
    cp true/true.rho inv/inv.rho \n
    cd .. 
    
    \b b) if necessary, run the script \b S0_prepare.py to update \b wavelet files, as
    
    # python /work/ess-leit/seisplatform/FD/SOFI2DAC/pyfiles/S0_prepare.py
    
    here the \b wavelet file was copied from dir. \b source.
    
    \b c) run the script \b cmds_true.sh to generate the data, as
    
    # bsub <cmds_true_p1.lsf
    # bsub <cmds_true_p2.lsf
    # bsub <cmds_true_p*.lsf
    
    
    @return observed data (fdio/sz***/*.su)
    
    here sz*** denotes the eventname, *.su the pressure data.
    """
    pass

def Sd_run_FWI_one_stage():
    """! Describe how to sumbit FWI jobs for one stage.
    In order to apply \b Multi-scale strategy for FWI, we filter 
    the seismograms with a lowpass window, such as Butterworth type. the
    corresponding filter parameters includes \b FC, \b ORDER, and 
    \b ZERO_PHASE which is defined in S0_prepare.py, the whole steps are:
    
    \b a) copy your \b current models into model/inv/, which is as
    
    cd model/ \n
    cp init/init.vp inv/inv.vp \n
    cp init/init.vs inv/inv.vs \n
    cp init/init.rho inv/inv.rho \n
    cd .. 
    
    \b b) if necessary, run the script \b S0_prepare.py to generate 
    all required commands, as
    
    # python /work/ess-leit/seisplatform/FD/SOFI2DAC/pyfiles/S0_prepare.py
    
    \b c) run the script \b cmds_filt.sh to filter the data and wavelets, as
    
    #only apply filtering at once.\n
    # python /work/ess-leit/seisplatform/FD/SOFI2DAC/pyfiles/Sb1_filt_allobs.py \n
    # python /work/ess-leit/seisplatform/FD/SOFI2DAC/pyfiles/Sb2_filt_allwavelet.py
    
    \b d) run optimization initialization, by
    
    # python /work/ess-leit/seisplatform/mod_OPTIM/init_optim_FWIac.py
    
    Before running, you should set the corresponding parameters for choosen optimization.
    
    \b e) run FWI in one stage by submitting the script \b submit_template.bash, 
    
    # bash submit_template.bash

    whose content is as
    
    #/usr/bin/bash \n
    #set iteration range\n
    mpiexe='/share/intel/2018u4/compilers_and_libraries_2018.5.274/linux/mpi/intel64/bin/mpiexec' \n
    spfroot='/work/ess-leit/seisplatform' \n
    iter_start=61 \n
    niter_max=90 \n
    nFD=8 \n
    ntest=3 \n
    #python ${spfroot}/mod_OPTIM/init_optim_FWIac.py \n
    mkdir -p log \n
     \n
    for ((iter=${iter_start}; iter<=${niter_max}; iter++)) \n
    do \n
        echo "at iteration $iter of ${niter_max}" \n
        #FDM\n
        for ((i=1;i<=${nFD};i++)) do \n
            bsub -J "GRAD${i}" <cmds_misf_grad_p${i}.lsf -K -o log/out.misf_grad_p${i}.it$iter -e log/err.misf_grad_p${i}.it$iter& \n
        done \n
        while [[ `ls log/out.misf_grad_p*.it$iter | wc -l` -ne ${nFD} ]] ; do \n
            echo "*** waiting GRAD ***" \n
            sleep 30 \n
        done \n
         \n
        id=$( expr ${nFD} + 1 ) \n
        bsub -J "GRAD" <cmds_misf_grad_p${id}.lsf -K -o log/out.misf_grad_p${id}.it$iter -e log/err.misf_grad_p${id}.it$iter \n
        while [[ `ls log/out.misf_grad_p${id}.it$iter | wc -l` -ne 1 ]] ; do \n
            echo "*** waiting GRAD2 ***" \n
            sleep 30 \n
        done \n
        #choose the descent \n
        #PNLCG \n
        python3 ${spfroot}/mod_OPTIM/PNLCG/init_PNLCG.py $iter \n
        python3 ${spfroot}/mod_OPTIM/PNLCG/descent_PNLCG.py \n
        #step length 0 \n
        python ${spfroot}/FD/SOFI2DAC/pyfiles/Sb4_sum_up_misfit.py test \n
        cp fdio/comm/misfit.dat optim/test/misfit0.dat \n
        python3 ${spfroot}/mod_OPTIM/linesearch/S0_alpha.py \n
        cp optim/test/m1.bin model/inv/inv.vp \n
        #step length 1\n
        for ((i=1;i<=${ntest};i++)) do \n
            bsub -J "test1-${i}" <cmds_misf_test_p${i}.lsf -K -o log/out.misf_test1_p${i}.it$iter -e log/err.misf_test1_p${i}.it$iter& \n
        done \n
        while [[ `ls log/out.misf_test1_p*.it$iter | wc -l` -ne ${ntest} ]] ; do \n
            echo "*** waiting test1 ***" \n
            sleep 30 \n
        done \n
        python ${spfroot}/FD/SOFI2DAC/pyfiles/Sb4_sum_up_misfit.py test \n
        cp fdio/comm/misfit.dat optim/test/misfit1.dat \n
        python3 ${spfroot}/mod_OPTIM/linesearch/S1_alpha.py \n
        cp optim/test/m2.bin model/inv/inv.vp \n
        \n
        #step length 2 \n
        for ((i=1;i<=${ntest};i++)) do \n
            bsub -J "test2-${i}" <cmds_misf_test_p${i}.lsf -K -o log/out.misf_test2_p${i}.it$iter -e log/err.misf_test2_p${i}.it$iter& \n
        done \n
        while [[ `ls log/out.misf_test2_p*.it$iter | wc -l` -ne ${ntest} ]] ; do \n
            echo "*** waiting test2 ***" \n
            sleep 30 \n
        done \n
        python ${spfroot}/FD/SOFI2DAC/pyfiles/Sb4_sum_up_misfit.py test \n
        cp fdio/comm/misfit.dat optim/test/misfit2.dat \n
        python3 ${spfroot}/mod_OPTIM/linesearch/S2_alpha.py \n
        cp optim/test/mod_optim.bin model/inv/inv.vp \n
         \n
    done \n
    @return inversion results (optim/iter*/mod.bin)
    
    here * denotes the iteration index.
    """
    pass

if __name__ == "__main__":
    Areadme()