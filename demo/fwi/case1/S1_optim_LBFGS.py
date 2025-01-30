#!/usr/bin/env python3
"""
the inversion example for Rosenbrock problem
Usage:

"""
import os
import sys
import glob
import numpy  as np
import pandas as pd

def split_jobs(jobs, cell, ncores):
  """
  split jobs into subjobs for one lsf file
  """
  nsrc = len(jobs)
  num_in_group = ncores//cell
  subjobs = [jobs[i:i+num_in_group] for i in range(0, nsrc, num_in_group)]
  return subjobs

def generate_cmdfile(jobnm, cmdlist, queuenm, ncores, ncore, wtime, headers, parallel=1):
  """
  generate jobs with lsf system
  if necessary, change corresponding parameters.
  """
  hdconf ='''#!/bin/bash
    #BSUB -q queuenm
    #BSUB -n ncores
    #BSUB -R "span[ptile=ncore]"
    #BSUB -W wtime
    #BSUB -e log/%J.err
    #BSUB -o log/%J.out

    module load intel/2020.4
    module load mpi/intel/2020.4
    module load fftw/3.3.8
    module load matlab/2017b
    
    '''
  hdconf = hdconf.replace('queuenm',queuenm)
  hdconf = hdconf.replace('ncores','%d'%ncores)
  hdconf = hdconf.replace('ncore','%d'%ncore)
  hdconf = hdconf.replace('wtime',wtime)
  hdconf = hdconf.replace("    ","")
    
  fp = open(jobnm,'w')
  # bash headers
  a = fp.write(hdconf)
  # environment variables
  for key, value in headers.items():
    s = key + '=' + str(value)
    a = fp.write("%s\n"%s)
  a = fp.write("\n")
  # write users' commands
  num = len(cmdlist)
  if parallel == 1:
    a = fp.write('commands=( \n')
    for i in range(num):
      a = fp.write('    "'+cmdlist[i] + ' &"\n')
    a = fp.write(') \n\n')
    # Run each command
    a = fp.write('# Record all job PIDs \n')
    a = fp.write('job_pids=() \n')
    a = fp.write('job_status=() \n')
    a = fp.write('for cmd in "${commands[@]}"; do\n')
    a = fp.write('    eval $cmd \n')
    a = fp.write('    job_pids+=($!) \n')
    a = fp.write('done\n\n')
    a = fp.write('# Use the wait command to wait for and check the status of each job \n')
    a = fp.write('for pid in "${job_pids[@]}"; do\n')
    a = fp.write('    wait $pid \n')
    a = fp.write('    job_status+=($?) \n')
    a = fp.write('done\n\n')
    # Check job statuses and print messages
    a = fp.write('# Check job statuses and print messages \n')
    a = fp.write('all_success=true \n')
    a = fp.write('for status in "${job_status[@]}"; do\n')
    a = fp.write('    if [ $status -ne 0 ]; then \n')
    a = fp.write('        all_success=false \n')
    a = fp.write('        break \n')
    a = fp.write('    fi \n')
    a = fp.write('done\n\n')
    # Print messages
    a = fp.write('if $all_success; then \n')
    a = fp.write('    echo "Exiting successfully." \n')
    a = fp.write('else \n')
    a = fp.write('    echo "Exiting unsuccessfully." \n')
    a = fp.write('fi \n\n')
  else:
    a = fp.write('commands=( \n')
    for i in range(num):
      a = fp.write('    "'+cmdlist[i] + '"\n')
    a = fp.write(') \n\n')
    # Run and check each Python command
    a = fp.write('for cmd in "${commands[@]}"; do\n')
    a = fp.write('    eval $cmd \n')
    a = fp.write('    if [ $? -ne 0 ]; then\n')
    a = fp.write('        echo "Exiting unsuccessfully." \n')
    a = fp.write('        echo "$cmd failed. Exiting..." >&2 \n')
    a = fp.write('        exit 1\n')
    a = fp.write('    fi\n')
    a = fp.write('done\n\n')
    a = fp.write('echo "Exiting successfully."\n')
  # close the file
  fp.close()


#INPut paras
fdir        =sys.argv[1] # optim
iter_start  =int(sys.argv[2]) # inversion starts from iter_start for this submission
iter_end    =int(sys.argv[3]) # inversion ends to iter_end for this submission
foptim      =os.path.join(fdir,'optim.json')
optim       =eval(open(foptim).read())
mpiexep     =optim['mpiexep']
mpiexec     =optim['mpiexec']
prjroot     =optim['prjroot']
optimroot   =optim['optimroot']
fwiroot     =optim['fwiroot']
fdroot      =optim['fdroot']
fevn_info   =optim['fevn_info']
ncores      =optim['ncores']
nproc       =optim['nproc']

py3         ='python3'
dirjob      ='job'
fdexe       =os.path.join(fdroot, 'bin', 'sofi2Dac')
adjexe      =os.path.join(fdroot, 'bin', 'sofi2Dac_grad')
npr         =nproc[0] * nproc[1]
ncore       =40 #ptile=ncore
queue_we      ='short'
queue_we2     ='ser'
queue_single  ='smp'
wtime         ='10:00'

headers = {
  "prjroot"   :prjroot,
  "fwiroot"   :fwiroot,
  "optimroot" :optimroot,
  "fdexe"     :fdexe,
  "adjexe"    :adjexe,
  "mpiexep"   :mpiexep,
  "mpiexec"   :mpiexec,
  "py3"       :py3
}

#source list
evn_info = pd.read_csv(fevn_info)
shot_list = list(evn_info['station'])
test = np.array(evn_info['test'])
shot_test_list = list(evn_info['station'][test==1])

"""
Initialization for current stage
"""
# make directory
for fnm in glob.glob('job/*'):
  os.remove(fnm)
for iterc in range(iter_start, iter_end+1, 1):
    os.makedirs(os.path.join(fdir, 'iter'+str(iterc)), exist_ok=True)

"""
Inversion loop
"""
df = pd.DataFrame(columns=['iter', 'step', 'inner', 'status', 'cmd', 'outfile', 'errfile'])
idx = -1
for iterc in range(iter_start, iter_end+1, 1):
  # print information
  print('*'*80)
  print('starts to generate the commands at iter. %4d of [%4d %4d]'%(iterc, iter_start, iter_end))  
  
  # run FD simulation
  jobs = []
  for evt in shot_list:
      cmd = 'cd '+os.path.join('${prjroot}', 'fdio', evt) + ' && ${mpiexec} -np '+ str(npr) + ' ${fdexe} grad.json'
      jobs.append(cmd)
  subjobs = split_jobs(jobs, npr, ncores)
  for i in range(len(subjobs)):
    jobnm = os.path.join(dirjob, 'FDa_p%02d.lsf'%(i+1))
    generate_cmdfile(jobnm, subjobs[i], queue_we, ncores, ncore, wtime, headers, 1)
    name = jobnm.split('/')[-1][:-4]
    out = 'log/out.%s.it%d'%(name, iterc)
    err = 'log/err.%s.it%d'%(name, iterc)
    cmd = 'bsub -J %s <%s -o %s -e %s >>log/lsf.log'%(name, jobnm, out, err)
    idx += 1
    df.loc[idx] = [iterc, 1, 0, 0, cmd, out, err]
  
  # collect the seismograms
  subjob = [
      '${mpiexep} -np 40 ${py3} ' + os.path.join('${fwiroot}', 'prob', 'S1_collect_seis.py')+' '+fdir+' all',
      ]
  jobnm = os.path.join(dirjob, 'FDb.lsf')
  generate_cmdfile(jobnm, subjob, queue_we, ncores, ncore, wtime, headers, 0)
  name = jobnm.split('/')[-1][:-4]
  out = 'log/out.%s.it%d'%(name, iterc)
  err = 'log/err.%s.it%d'%(name, iterc)
  cmd = 'bsub -J %s <%s -o %s -e %s >>log/lsf.log'%(name, jobnm, out, err)
  idx += 1
  df.loc[idx] = [iterc, 2, 0, 0, cmd, out, err]
  
  # calculate adjoint source
  subjob = [
      '${mpiexep} -np 40 ${py3} ' + os.path.join('${fwiroot}', 'prob', 'S2_process_seis.py')+' '+fdir,
      '${mpiexep} -np 40 ${py3} ' + os.path.join('${fwiroot}', 'prob', 'S3_adjoint_src.py')+' '+fdir+' 1'
      ]
  jobnm = os.path.join(dirjob, 'adj_src.lsf')
  generate_cmdfile(jobnm, subjob, queue_we, ncores, ncore, wtime, headers, 0)
  name = jobnm.split('/')[-1][:-4]
  out = 'log/out.%s.it%d'%(name, iterc)
  err = 'log/err.%s.it%d'%(name, iterc)
  cmd = 'bsub -J %s <%s -o %s -e %s >>log/lsf.log'%(name, jobnm, out, err)
  idx += 1
  df.loc[idx] = [iterc, 3, 0, 0, cmd, out, err]
  
  # run BD simulation
  jobs = []
  for evt in shot_list:
      cmd = 'cd '+os.path.join('${prjroot}', 'fdio', evt) + ' && ${mpiexec} -np '+ str(npr) + ' ${adjexe} grad.json'
      jobs.append(cmd)
  subjobs = split_jobs(jobs, npr, ncores)
  for i in range(len(subjobs)):
    jobnm = os.path.join(dirjob, 'BD_p%02d.lsf'%(i+1))
    generate_cmdfile(jobnm, subjobs[i], queue_we, ncores, ncore, wtime, headers, 1)
    name = jobnm.split('/')[-1][:-4]
    out = 'log/out.%s.it%d'%(name, iterc)
    err = 'log/err.%s.it%d'%(name, iterc)
    cmd = 'bsub -J %s <%s -o %s -e %s >>log/lsf.log'%(name, jobnm, out, err)
    idx += 1
    df.loc[idx] = [iterc, 4, 0, 0, cmd, out, err]
  
  # calculate gradient
  subjob = [
      '${py3} ' + os.path.join('${fwiroot}', 'prob', 'S4_calc_gradient.py')+' '+fdir,
      '${py3} ' + os.path.join('${fwiroot}', 'prob', 'S5_process_gradient.py')+' '+fdir
      ]
  jobnm = os.path.join(dirjob, 'gradient.lsf')
  generate_cmdfile(jobnm, subjob, queue_we, ncores, ncore, wtime, headers, 0)
  name = jobnm.split('/')[-1][:-4]
  out = 'log/out.%s.it%d'%(name, iterc)
  err = 'log/err.%s.it%d'%(name, iterc)
  cmd = 'bsub -J %s <%s -o %s -e %s >>log/lsf.log'%(name, jobnm, out, err)
  idx += 1
  df.loc[idx] = [iterc, 5, 0, 0, cmd, out, err]
  
  # calculate the descent direction
  subjob = [
      '${py3} ' + os.path.join('${optimroot}', 'LBFGS', 'init.py')+' '+fdir+' '+str(iterc), 
      '${py3} ' + os.path.join('${optimroot}', 'LBFGS', 'descent.py')+' '+fdir
      ]
  jobnm = os.path.join(dirjob, 'descent_'+str(iterc)+'.lsf')
  generate_cmdfile(jobnm, subjob, queue_single, 1, 1, wtime, headers, 0)
  name = jobnm.split('/')[-1][:-4]
  out = 'log/out.%s.it%d'%(name, iterc)
  err = 'log/err.%s.it%d'%(name, iterc)
  cmd = 'bsub -J %s <%s -o %s -e %s >>log/lsf.log'%(name, jobnm, out, err)
  idx += 1
  df.loc[idx] = [iterc, 6, 0, 0, cmd, out, err]
  
  # calculate the step length
  ## collect current misfit
  subjob = [
      '${mpiexep} -np 40 ${py3} ' + os.path.join('${fwiroot}', 'prob', 'S3_adjoint_src.py')+' '+fdir+' 3',
      '${py3} '+ os.path.join('${optimroot}', 'linesearch', 'S0_alpha.py')+' '+fdir
      ]
  jobnm = os.path.join(dirjob, 'misfit0.lsf')
  generate_cmdfile(jobnm, subjob, queue_we, ncores, ncore, wtime, headers, 0)
  name = jobnm.split('/')[-1][:-4]
  out = 'log/out.%s.it%d'%(name, iterc)
  err = 'log/err.%s.it%d'%(name, iterc)
  cmd = 'bsub -J %s <%s -o %s -e %s >>log/lsf.log'%(name, jobnm, out, err)
  idx += 1
  df.loc[idx] = [iterc, 7, 0, 0, cmd, out, err]
  
  ## calculate misfit at 1st test
  ### FD
  jobs = []
  for evt in shot_test_list:
      cmd = 'cd '+os.path.join('${prjroot}', 'fdio', evt) + ' && ${mpiexec} -np '+ str(npr) + ' ${fdexe} test.json'
      jobs.append(cmd)
  subjobs = split_jobs(jobs, npr, ncores)
  for i in range(len(subjobs)):
    jobnm = os.path.join(dirjob, 'misfit1a_p%02d.lsf'%(i+1))
    generate_cmdfile(jobnm, subjobs[i], queue_we, ncores, ncore, wtime, headers, 1)
    name = jobnm.split('/')[-1][:-4]
    out = 'log/out.%s.it%d'%(name, iterc)
    err = 'log/err.%s.it%d'%(name, iterc)
    cmd = 'bsub -J %s <%s -o %s -e %s >>log/lsf.log'%(name, jobnm, out, err)
    idx += 1
    df.loc[idx] = [iterc, 8, 0, 0, cmd, out, err]
  
  ### collect misfit
  subjob = [
      '${mpiexep} -np 40 ${py3} ' + os.path.join('${fwiroot}', 'prob', 'S1_collect_seis.py')+' '+fdir+' test',
      '${mpiexep} -np 40 ${py3} ' + os.path.join('${fwiroot}', 'prob', 'S2_process_seis.py')+' '+fdir,
      '${mpiexep} -np 40 ${py3} ' + os.path.join('${fwiroot}', 'prob', 'S3_adjoint_src.py')+' '+fdir+' 2',
      '${py3} '+ os.path.join('${optimroot}', 'linesearch', 'S1_alpha.py')+' '+fdir
      ]
  jobnm = os.path.join(dirjob, 'misfit1b.lsf')
  generate_cmdfile(jobnm, subjob, queue_we, ncores, ncore, wtime, headers, 0)
  name = jobnm.split('/')[-1][:-4]
  out = 'log/out.%s.it%d'%(name, iterc)
  err = 'log/err.%s.it%d'%(name, iterc)
  cmd = 'bsub -J %s <%s -o %s -e %s >>log/lsf.log'%(name, jobnm, out, err)
  idx += 1
  df.loc[idx] = [iterc, 9, 0, 0, cmd, out, err]
  
  ## calculate misfit at 2nd test
  ### FD
  jobs = []
  for evt in shot_test_list:
      cmd = 'cd '+os.path.join('${prjroot}', 'fdio', evt) + ' && ${mpiexec} -np '+ str(npr) + ' ${fdexe} test.json'
      jobs.append(cmd)
  subjobs = split_jobs(jobs, npr, ncores)
  for i in range(len(subjobs)):
    jobnm = os.path.join(dirjob, 'misfit2a_p%02d.lsf'%(i+1))
    generate_cmdfile(jobnm, subjobs[i], queue_we, ncores, ncore, wtime, headers, 1)
    name = jobnm.split('/')[-1][:-4]
    out = 'log/out.%s.it%d'%(name, iterc)
    err = 'log/err.%s.it%d'%(name, iterc)
    cmd = 'bsub -J %s <%s -o %s -e %s >>log/lsf.log'%(name, jobnm, out, err)
    idx += 1
    df.loc[idx] = [iterc, 10, 0, 0, cmd, out, err]
  
  ### collect misfit
  subjob = [
      '${mpiexep} -np 40 ${py3} ' + os.path.join('${fwiroot}', 'prob', 'S1_collect_seis.py')+' '+fdir+' test',
      '${mpiexep} -np 40 ${py3} ' + os.path.join('${fwiroot}', 'prob', 'S2_process_seis.py')+' '+fdir,
      '${mpiexep} -np 40 ${py3} ' + os.path.join('${fwiroot}', 'prob', 'S3_adjoint_src.py')+' '+fdir+' 2',
      '${py3} '+ os.path.join('${optimroot}', 'linesearch', 'S2_alpha.py')+' '+fdir
      ]
  jobnm = os.path.join(dirjob, 'misfit2b.lsf')
  generate_cmdfile(jobnm, subjob, queue_we, ncores, ncore, wtime, headers, 0)
  name = jobnm.split('/')[-1][:-4]
  out = 'log/out.%s.it%d'%(name, iterc)
  err = 'log/err.%s.it%d'%(name, iterc)
  cmd = 'bsub -J %s <%s -o %s -e %s >>log/lsf.log'%(name, jobnm, out, err)
  idx += 1
  df.loc[idx] = [iterc, 11, 0, 0, cmd, out, err]
  
  # save current information
  subjob = [
      '${py3} ' + os.path.join('${fwiroot}', 'prob', 'S6_copy_files.py')+' '+fdir,
      '${py3} ' + os.path.join('${fwiroot}', 'prob', 'S7_check_stopping.py')+' '+fdir
      ]
  jobnm = os.path.join(dirjob, 'save.lsf')
  generate_cmdfile(jobnm, subjob, queue_single, 1, 1, wtime, headers, 0)
  name = jobnm.split('/')[-1][:-4]
  out = 'log/out.%s.it%d'%(name, iterc)
  err = 'log/err.%s.it%d'%(name, iterc)
  cmd = 'bsub -J %s <%s -o %s -e %s >>log/lsf.log'%(name, jobnm, out, err)
  idx += 1
  df.loc[idx] = [iterc, 12, 0, 0, cmd, out, err]
  print('finish optimization at iter. %4d of [%4d %4d] \n'%(iterc, iter_start, iter_end))

#write the tasks
jobNM = os.path.join(dirjob, 'job.csv')
df.to_csv(jobNM, index=False)
  