# Full waveform inversion Documentation

>Full Waveform Inversion (FWI) is a sophisticated seismic imaging technique that achieves highly accurate subsurface models by leveraging complete waveform data. This comprehensive approach results in detailed and precise subsurface images, providing high-resolution imaging of structures crucial for applications such as oil and gas exploration, and environmental studies. FWI accurately captures and reproduces the complex wavefield interactions within the Earth's subsurface. However, this process demands significant computational power due to the intricate wave equations and the vast amounts of data involved. Achieving accurate results typically involves numerous iterations, making the overall inversion process time-consuming and computationally intensive. Our optimization software includes a two-dimensional acoustic Full Waveform Inversion (FWI) implementation to demonstrate its applicability to solving complex geophysical problems. The exceptional results suggest that our optimization framework is highly effective in addressing intricate inverse challenges. We use the LSF management system to manage a large number of computational tasks. If a job is interrupted unfortunately, our inversion process can detect the point at which the previous task was terminated and will continue to complete the remaining tasks. This also highlights the advantage of our software's reliance on file-based data transfer. This documentation aims to guide you on how to organize the optimization workflow for FWI study and submit jobs to an HPC cluster. The optimization toolbox [optim](https://github.com/beyondTaoLei/optim) has been developed by Tao Lei (leit@sustech.edu.cn) and Wei Zhang (zhangwei@sustech.edu.cn). If you encounter any issues while using this software, please do not hesitate to contact us.
***

## 0. Demo structure
The root directory contains several key components essential for understanding and utilizing the provided resources. It includes this documentation, the workflow example (**case1**), the collection of main functions (**prob**), some drawing scripts located in **plot**, modified seismic forward & adjoint operators **SOFI2DAC** and the associated template input file in **conf** folder. 

Most parameters will be set using the script **S0_prepare.py** within the workflow example(**case1**). The following directory tree presents the main structure for solving 2D acoustic full waveform inversion problem:

```
|-- README.md                   >>> this documentation
|-- case1                       >>> workflow example  
|   |-- P0_prepare_files.py      
|   |-- P1_prepare_json.py
|   |-- P2_prepare_wavelet.py
|   |-- P3_generate_data.py
|   |-- S0_prepare.py           >>> prepare the global JSON file for 1st stage
|   |-- S0_prepare_fc5.py       >>> ... for 2nd stage
|   |-- S0_prepare_fc7.py       >>> ... for 3rd stage
|   |-- S1_optim_LBFGS.py
|   |-- S2_run.py
|   `-- model                   >>> initial and true Marmousi II models
|-- plot                        >>> figure presentation
|   |-- plot_misfit.py
|   |-- plot_mod.py
|   `-- plot_seis.py
|-- prob                        >>> main functions used in workflow
|   |-- S1_collect_seis.py
|   |-- S2_process_seis.py
|   |-- S3_adjoint_src.py
|   |-- S4_calc_gradient.py
|   |-- S5_process_gradient.py
|   |-- S6_copy_files.py
|   |-- S7_check_stopping.py
|   |-- misfit_funs.py
|   |-- update_optim.py
|   |-- util.py
|   `-- wavelet.py
|-- SOFI2DAC                    >>> external seismic forward & adjoint operators
|   |-- bin
|   |-- contrib
|   |-- mfiles
|   `-- src
`-- conf                        
    `-- FD.json                 >>> template input file for SOFI2DAC software

      Figure 1. The structure of full waveform inversion demo
```
## 1. Configure environment

Our software encompasses both parallel programs in Python and C. Given that the environment management tool, *conda*, often prioritizes the Python MPI parallel launcher in Linux system environment variables, we recommend compiling the seismic programs written in C before activating the conda virtual environment.

### 1.1 Compile seismic modeling operators
We first utilize the [sofi2d](https://gitlab.kit.edu/kit/gpi/ag/software/sofi2d) software for acoustic wave-equation simulation, and then apply the adjoint state method to calculate the gradient for FWI study. These two operators in our software are developed based on the original [sofi2d](https://gitlab.kit.edu/kit/gpi/ag/software/sofi2d) software package. For detailed instructions on compiling and usage, please refer to the *sofi2d* documentation. Now, let's outline the compilation process for our seismic modeling operators.
- Ensure that an MPI compiler for C, *mpicc*, is installed, and compile the library for the input and output of seismic modeling operators
  ```sh
  # adjust your path
  cd /path/to/your/optim/software/demo/fwi/SOFI2DAC/src/IN_OUT
  make lib
  ```

- Compile the 2D acoustic wave-equation finite-difference simulation operator and its gradient operator
  ```sh
  # change the directory
  cd ..
  # generate forward modeling operator
  make sofi2Dac
  # generate gradient operator
  make sofi2Dac_grad
  # check the result of compilation
  cd ../bin && ls -rlt
  ```
Upon completion of the compilation, you will find two operators: *sofi2Dac* and *sofi2Dac_grad*. These operators will be utilized within our FWI workflow.

### 1.2 Configure Python environment
To ensure your software runs smoothly, install the necessary Python libraries by following these steps:
- Create a new environment (optional): It's recommended to create and use a virtual environment to manage dependencies. Replace *myenv* with your desired environment name.
  ```sh
  conda create -n myenv
  conda activate myenv
  ```
- Install dependencies: Use the following command to install the required libraries.
  ```sh
  conda install -c conda-forge numpy scipy pandas matplotlib mpi4py obspy
  ```
- List installed libraries: To verify that the libraries are installed successfully, list the installed packages.
  ```sh
  conda list
  ```
By following these steps, you'll ensure that all necessary packages are installed and your Python environment is set up properly.

## 2. Prepare observed data

- Specify the path to the FWI demo and name your project in the terminal:
  ```sh
  demoroot=/path/to/software/optim/demo/fwi
  mkdir prj
  cd prj
  ```
  Here, **demoroot** indicates the root path of the example codes. Specifying this path in the terminal will facilitate copying and using the subsequent codes for study.

- Modify the program path variables in the script **${demoroot}/case1/S0_prepare.py**  
  ```sh
  'Program path'      :'comment',
  'optimroot'         :'/path/to/software/optim',
  'mpiexec'           :'/path/to/your/C/MPI/parallel/launcher',       # c mpirun
  'mpiexep'           :'/path/to/your/python/MPI/parallel/launcher',  # python mpirun
  ```
  **Parameters**
  - ***optimroot: string***  
  The root path to the optimization toolbox codes. 
  - ***mpiexec: string***  
  The path to C MPI parallel launcher. 
  - ***mpiexep: string***  
  The path to Python MPI parallel launcher.

- Write the global parameters for the whole project
  ```sh
  python ${demoroot}/case1/S0_prepare.py optim
  ```
  :bell: Output: **optim/optim.json**
  We place all optimization parameters in the script **S0_prepare.py**. Each program in the inverse workflow will extract the necessary parameters from the output JSON file. The input parameter **optim** in the command specifies the output folder to store the global parameters.

- Copy the Marmousi models from the *optim* package into the current project directory, distribute the true model to a specific location, and display the model as
  ```sh
  cp -r ${demoroot}/case1/model .
  cp model/true/true.vp model/inv/inv.vp
  cp model/true/true.rho model/inv/inv.rho
  python ${demoroot}/plot/plot_mod.py optim model/inv/inv.vp
  ```
  | Target          | :bell: Output |
  | ------          | ------ |
  | vp model        | **model/inv/inv.vp** |
  | density model   | **model/inv/inv.rho** |

- Prepare the geometry for the current study, which includes 80 sources and 801 receivers located on the surface. 
  ```sh
  python ${demoroot}/case1/P0_prepare_files.py optim
  ``` 
  | Target              | :bell: Output |
  | ------              | ------ |
  | source information  | **list/eventinfo.csv** |
  | source wavelet      | **source/srcname.su** |
  | receiver position   | **fdio/srcname/receiver.dat** |
  | gradient taper      | **list/taper.bin** |

  where the keyword **srcname** indicates the source name at each source position. the detailed source information is provided by the file **eventinfo.csv**, including the name, position, and delayed time, etc. 

- Prepare the input json file of *SOFI2DAC* software for each source
    ```sh
    python ${demoroot}/case1/P1_prepare_json.py optim
    ```
    | Target                  | :bell: Output |
    | ------                  | ------ |
    | forward modeling        | **fdio/srcname/test.json** |
    | gradient calculation    | **fdio/srcname/grad.json** |

- Write the excited source wavelet signal for each source
    ```sh
    python ${demoroot}/case1/P2_prepare_wavelet.py optim 0
    ```
    | Target            | :bell: Output |
    | ------            | ------ |
    | source signal     | **fdio/srcname/wavelet.su** |

- Generate the observed data for all shots
  ```sh
  # first generate the commands
  python ${demoroot}/case1/P3_generate_data.py optim
  # and then execute with platform LSF
  python ${demoroot}/case1/S2_run.py optim
  ```
  | Target              | :bell: Output |
  | ------              | ------ |
  | job list            | **job/job.csv** |
  | observed data       | **fdio/srcname/seis_p_shot1.su** |
  
  where the first script will generate the forward modeling commands and store them in the file **job/job.csv**, while the second script will execute these commands and generate the observed seismograms.
  **Note that**, if you do not wish to use the LSF job management system, you should modify the function *generate_cmdfile* in the script *case1/P3_generate_data.py*.

- Display the observed data
  ```sh
  python ${demoroot}/plot/plot_seis.py data/obs0.h5 sz001.P
  ```

## 3. Run nonlinear inverse workflow

- Prepare the initial model to start the inversion process
  ```sh
  cp model/init/init.vp model/inv/inv.vp
  cp model/init/init.rho model/inv/inv.rho
  python ${demoroot}/plot/plot_mod.py optim model/inv/inv.vp
  ```
  | Target          | :bell: Output |
  | ------          | ------ |
  | vp model        | **model/inv/inv.vp** |
  | density model   | **model/inv/inv.rho** |

- Write the global parameters for the whole project
  ```sh
  python ${demoroot}/case1/S0_prepare.py LBFGS
  ```
  :bell: Output: **LBFGS/optim.json**
  We place all optimization parameters in the script **S0_prepare.py**. Each program in the inverse workflow will extract the necessary parameters from the output JSON file. The input parameter **LBFGS** in the command specifies the output folder to store the global parameters.

- Generate the LBFGS workflow from iteration 1 to 30
  ```sh
  python ${demoroot}/case1/S1_optim_LBFGS.py LBFGS 1 30
  ```
  :bell: Output: **job/job.csv**
  where the file **job.csv** contains the commands for the inverse workflow.

- Implement data pre-processing for observed data, for example
  ```sh
  grep -v "adjoint_src" job/adj_src.lsf > job/process.lsf
  cp data/obs0.h5 data/syn0.h5
  bsub <job/process.lsf
  mv data/syn.h5 data/obs.h5
  ```
  :bell: Output: **data/obs.h5**

- Write the excited source wavelet signal for each source
  ```sh
  # filtering the wavelet rather than the synthetic data
  python ${demoroot}/case1/P2_prepare_wavelet.py LBFGS 1
  # update the optim.json without filtering the synthetic data
  python ${demoroot}/prob/update_optim.py LBFGS timefilt 0
  ```
  | Target                | :bell: Output |
  | ------                | ------ |
  | source signal         | **fdio/srcname/wavelet.su** |
  | global parameters     | **LBFGS/optim.json** |
  
  We choose to apply the lowpass filtering to the source wavelet instead of the synthetic data. The advantage of this approach is that it not only ensures records have consistent processing effects, but also prevents the forward wavefield from containing high-frequency components.

- Execute the LBFGS inverse workflow
  ```sh
  python ${demoroot}/case1/S2_run.py LBFGS
  ```
  | Target              | :bell: Output |
  | ------              | ------ |
  | job list            | **job/job.csv** |
  | evaluated model     | **LBFGS/iter#/mod.bin** |
  | data misfit         | **LBFGS/iter#/fcost.dat** |
  | gradient            | **LBFGS/iter#/grad.bin** |

  The file **job.csv** contains the commands for the inverse workflow, and it updates after finishing each command. There are three different statuses for each command: 0 means the task has not started, 1 means the task is finished, and -1 means the task failed. You need to find the corresponding output file and error log to view the specific reason. **mod.bin** represents the evaluated model at the **#**-th iteration, and **fcost.dat** and **grad.bin** represents the total data misfit and the model gradient.

- Repeat the above steps for the following stages
  ```sh
  # 2nd stage at fc = 5.0 Hz
  ## write the global parameters (specially)
  python ${demoroot}/case1/S0_prepare_fc5.py LBFGS
  ## generate the LBFGS workflow (specially)
  python ${demoroot}/case1/S1_optim_LBFGS.py LBFGS 31 60
  ## implement data pre-processing for observed data
  grep -v "adjoint_src" job/adj_src.lsf > job/process.lsf
  cp data/obs0.h5 data/syn0.h5
  bsub <job/process.lsf
  mv data/syn.h5 data/obs.h5
  ## filter the wavelet rather than the synthetic data
  python ${demoroot}/case1/P2_prepare_wavelet.py LBFGS 1
  ## update the optim.json without filtering the synthetic data
  python ${demoroot}/prob/update_optim.py LBFGS timefilt 0
  ## execute the workflow
  python ${demoroot}/case1/S2_run.py LBFGS

  # 3rd stage at fc = 7.0 Hz
  ## write the global parameters (specially)
  python ${demoroot}/case1/S0_prepare_fc7.py LBFGS
  ## generate the LBFGS workflow (specially)
  python ${demoroot}/case1/S1_optim_LBFGS.py LBFGS 61 90
  ## implement data pre-processing for observed data
  grep -v "adjoint_src" job/adj_src.lsf > job/process.lsf
  cp data/obs0.h5 data/syn0.h5
  bsub <job/process.lsf
  mv data/syn.h5 data/obs.h5
  ## filter the wavelet rather than the synthetic data
  python ${demoroot}/case1/P2_prepare_wavelet.py LBFGS 1
  ## update the optim.json without filtering the synthetic data
  python ${demoroot}/prob/update_optim.py LBFGS timefilt 0
  ## execute the workflow
  python ${demoroot}/case1/S2_run.py LBFGS
  ```

- Display the inverse results using the LBFGS approach
  ```sh
  # plot model
  python ${demoroot}/plot/plot_mod.py LBFGS LBFGS/iter90/mod.bin
  # plot data
  python ${demoroot}/plot/plot_seis.py data/obs0.h5 sz001 P
  python ${demoroot}/plot/plot_seis.py LBFGS/iter90/syn.h5 sz001 P
  # plot misfit func.
  python ${demoroot}/plot/plot_misfit.py LBFGS
  ```

## 4. Set json Parameters for nonlinear inverse algorithm

In this section, we will discribe how to set the parameters in the configuration JSON file. The final file should include the following key components: 1) program path module, 2) file path information, 3) forward modeling module, 4) date manipulation module, 5) gradient calculation module, 6) descent direction module, 7) model update module, 8) inversion loop. We will go over each of these components in detail.

### 4.1 Program path module
  Since we run all commands from the project root path, we need to define the path to the optimization software, C MPI parallel launcher and Python MPI parallel launcher, as
  ```sh
  'Program path'      :'comment',
  'optimroot'         :'/path/to/software/optim',
  'mpiexec'           :'/path/to/your/C/MPI/parallel/launcher',       # c mpirun
  'mpiexep'           :'/path/to/your/python/MPI/parallel/launcher',  # python mpirun
  ```
  **Parameters**
  - ***optimroot: string***  
  The root path to the optimization toolbox codes. 
  - ***mpiexec: string***  
  The path to C MPI parallel launcher. 
  - ***mpiexep: string***  
  The path to Python MPI parallel launcher. 

### 4.2 File path information
  We build the inverse workflow with a series of procedures, including the geophysical modules and the optimization modules, and each procedure connects the predecessors and successors with files. Therefore, it is necessary to specify the locations of the model, cost function, gradient, and other relevant files. The contents of these files may be updated over iterations.
  ```sh
  'File path'         :'comment',
  'fevn_info'         :'list/eventinfo.csv', 
  'fmodfd'            :'model/inv/inv.vp',    
  'fmod'              :'model/inv/inv.vp',    # optim needs
  'fcost'             :'list/fcost.dat',      # optim needs
  'fgrad'             :'list/grad.bin',       # optim needs
  'grad_taper_file'   :'list/taper.bin',      
  ```
  **Parameters**
  - ***fevn_info: input, string***  
    The path for storing the source information, where the data is saved in *csv* format. The key 'station' in the header indicates the source name, while 'x' and 'z' indicate the coordinates of the source, 'tshift' indicate the time delay, and 'fc' denotes the center frequency of the wavelet.
  - ***fmodfd: input&output, string***  
    The path for storing the current vp model, where the data is saved in a 32-bit floating-point binary file. The input data is arranged in *column-major* order, with the coordinate origin at the top-left corner. The x-axis extends to the right, and the y-axis extends downward. Note that, the associated density model, named with the suffix .rho, has also to be provided.
  - ***fmod: output, string***  
    The path for storing the resampled output of the model is within the file ***fmodfd***. For this study, we do not resample the grids during calculating the gradient.
  - ***fcost: output, string*** 
    The path for storing the misfit function, where the data is saved in *ASCII* floating-point type.
  - ***fgrad: output, string***  
    The path for storing the gradient, where the data is saved in a 32-bit floating-point binary file.
  - ***grad_taper_file: input, string***  
    The path to user-defined taper file for the gradient, where the data is saved in a 32-bit floating-point binary file.

### 4.3 Forward modeling module
  Next, we will detail the parameters specifically for the forward modeling. 
  ```sh
  'Forward part'      :'comment',
  'ncores'            :40,                      
  'nproc'             :[4, 1],                 
  'nxzfd'             :[1700, 350],            
  'oxzfd'             :[0.0, 0.0],              
  'dxzfd'             :[10.0, 10.0],           
  'nt'                :6500,                    
  'delta'             :1.0e-3,                  
  'tshift_adj'        :0.2,                     
  'freesurface'       :1,                       
  'pml'               :15,                      
  'pmlvel'            :3000.0,                  
  'comp'              :['P'],                   
  ```
  **Parameters**
  - ***ncores: int***  
    The number of CPU cores within a single node. Since the tasks are submitted using the LSF management system, it is essential to declare this configuration.
  - ***nproc: list of int***  
    The number of MPI processors used for forward modeling in the x- and z-directions, respectively.
  - ***nxzfd: list of int***  
    The number of grid points of the input model for forward modeling in the x- and z-directions, respectively.
  - ***oxzfd: list of float***  
    The original coordinates of the input model for forward modeling in the x- and z-directions, respectively.
  - ***dxzfd: list of float***  
    The grid spacing of the input model for forward modeling in the x- and z-directions, respectively. For the *SOFI2DAC* software, the grid spacing is uniform in both directions.
  - ***nt: int***  
    The total number of time steps for the simulation.
  - ***delta: float***  
    The time stepping interval in second.
  - ***tshift_adj: float***  
    The deadline for backward modeling. Since there is a time delay when exciting the source, it is not necessary to calculate the correlation of the forward wavefield and the backward wavefield during this period. It is recommended to set this to the time delay of the source signal.
  - ***freesurface: int***  
    Whether to apply a plane stress free surface at the top of the global grid, as proposed by Levander (1988). A value of 1 means yes, while 0 means no.
  - ***pml: int***  
    The number of Perfectly Matched Layers (PML) applied as boundary conditions on each side face and the bottom face of the model grid. If the *freesurface* is set to 1, the boundary condition will also be applied on the top face. A width of the absorbing frame between 10-20 grid points should be sufficient.
  - ***pmlvel: float***  
    The attenuation velocity within the PML, measured in m/s. This value should be approximately equal to the propagation velocity of the dominant wave near the model boundaries.
  - ***comp: list of string***  
    The type of measurement components. For acoustic waveform inversion, only the pressure component ('P') is available.
  
### 4.4 Date manipulation module
  Next, we will detail the parameters for data processing. You can apply various date manipulation techniques to the seismograms, including time-domain integration, direct wave muting, filtering, trace killing, and time windowing.
  ```sh
  'Data processing'           :'comment',                 
  'velocity'                  :1,           
  'dirmute'                   :0,   
  'dirmute_taperlen'          :0.4,                     
  'dirmute_t0'                :0.20,                    
  'dirmute_lower'             :0.0, 
  'dirmute_upper'             :1700.0,                  
  'timefilt'                  :1,                       
  'fc'                        :2.5,                     
  'order'                     :6,                       
  'zero_phase'                :0,                       
  'trkill'                    :0,                       
  'trkill_lower'              :20.0,                    
  'trkill_upper'              :10000.0,                 
  'timewin'                   :1,                       
  'timewin_lower'             :0.1,                     
  'timewin_upper'             :10.0,          
  ```
  **Parameters**
  - ***velocity: int***  
    Whether to apply integration in the time domain. A value of 0 means integration is applied, while 1 means no integration.
  - ***dirmute: int***  
    Whether to apply direct wave muting. A value of 0 means no muting, 1 means direct waves will be muted according to the propagation speed of seismic waves, 2 means direct waves will be muted using *f-k* filtering, and 3 means direct waves will be muted using a user-defined muting file based on observed arrival travel times.
  - ***dirmute_taperlen: float***  
    The taper length when *dirmute* is set to 1. The muting time window starts from the direct waves and extends with a length defined by *dirmute_taperlen*.
  - ***dirmute_t0: float***  
    Time offset caused by the delayed excitation of the source when *dirmute* is set to 1.
  - ***dirmute_lower: float***  
    The propagating speed of the direct wave when *dirmute* is set to 1, or the stop band low corner speed when *dirmute* is set to 2.
  - ***dirmute_upper: float***  
    The stop band high corner speed when *dirmute* is set to 2. The data with propagating speeds between *dirmute_lower* and *dirmute_upper* will be muted.
  - ***timefilt: int***  
    Whether to apply lowpass filtering. A value of 1 means yes, while 0 means no.
  - ***fc: float***  
    The low cut-off frequency in Hz.
  - ***order: int***  
    The filter corners or order, with values of 4 or 6 being recommended based on experience.
  - ***zero_phase: int***  
    Whether to keep zero phase shift when filtering. A value of 1 means *true*, while 0 means *false*.
  - ***trkill: int***  
    Whether to apply trace killing. A value of 1 means yes, while 0 means no.
  - ***trkill_lower: float***  
    The minimum offset when *trkill* is set to 1. The traces outside the range [*trkill_lower*, *trkill_upper*] will be muted.
  - ***trkill_upper: float***  
    The maximum offset when *trkill* is set to 1.
  - ***timewin: int***  
    Whether to apply time windowing. A value of 1 means yes, while 0 means no.
  - ***timewin_lower: float***  
    The starting time of the window in seconds. Seismograms outside the time range [*timewin_lower*, *timewin_upper*] will be removed.
  - ***timewin_upper: float***  
    The end time of the window in seconds.
  
### 4.5 Gradient calculation module
  Now we will detail the parameters specifically for gradient calculations. 
  ```sh
  'Gradient part'             :'comment',                   
  'Lnorm'                     :2,                       
  'epsilon'                   :0.03,                    
  'grad_taper_shot'           :1,                       
  'srtradius1'                :1,                       
  'srtradius2'                :10,                      
  'grad_filt'                 :5,                       
  'grad_taper'                :1,                       
  'grad_taper_file'           :'list/taper.bin', # already shown in [File path information] module
  ```
  **Parameters**
  - ***Lnorm: int***  
    The type of misfit function. For now, you can choose 2 for the *L2* norm or 3 for *cross-correlation*.
  - ***epsilon: int***  
    The water level to stabilize the preconditioning operator, with a suggested range of 0 to 0.05.
  - ***grad_taper_shot: int***  
    Whether to apply a frustal taper around the source and receiver position for the gradient per shot. A value of 1 means yes, while 0 means no.
  - ***srtradius1: int***  
    The distance between the vertex of the taper and the source or receiver position when *grad_taper_shot* is set to 1. The frustal taper around the source/receiver positions decreases from a value of one at the edge of the taper (*r=srtradius2*) to a value of zero near the associated source/receiver positions (*r=srtradius1*).
  - ***srtradius2: int***   
    The distance between the edge of the taper and the source or receiver position when *grad_taper_shot* is set to 1.
  - ***grad_filt: int***  
    The smoothing size for the gradient in grid points. A value of 0 means no smoothing is applied.
  - ***grad_taper: int***   
    Whether to apply a taper to the final gradient. A value of 1 means yes, while 0 means no.
  - ***grad_taper_file: int***   
    The path to the user-defined taper file for the gradient when *grad_taper* is set to 1. The data should be saved in a 32-bit floating-point binary file.

### 4.6 Descent direction module
  For most optimization algorithms, users need to input few additional parameters to constrain the calculation of the descent direction. Here we will introduce the basic parameters for *l-BFGS* optimization method in the *OPTIM* software. In general, the *l*-BFGS method utilizes the gradients and models of the most recent iterations to construct the descent direction, while discarding information from earlier iterations to save storage.
  ```sh
  'l-BFGS'            :'comment',
  'l'                 :5          # optim needs
  ```
  **Parameters**
  - ***l: int***  
    The number of the most recent correction pairs $\{s_k, y_k\}$. Users can specify the number *l* according to the degree of complexity of the inverse problem. It is suggested to be less than 10.
  
### 4.7 Model update module
  Once the descent direction $p_k$ is generated by the above nonlinear optimization algorithms, line search strategy can help to find the optimal step length $\alpha_k$ for the model updating. We choose the parabolic fitting method to achieve this goal, and the relevant parameters are as follows:
  ```sh
  'Steplength'        :'comment',
  'try_old_sl'        :1,         # optim needs                
  'max_m0'            :4700.0,    # optim needs             
  'eps_scale'         :0.002,     # optim needs             
  'alpha_lb'          :0.1,       # optim needs             
  'alpha_ub'          :10,        # optim needs            
  'mod_limits'        :1,         # optim needs            
  'mod_lb'            :1000,      # optim needs            
  'mod_ub'            :5000,      # optim needs           
  ```
  **Parameters**
  - ***try_old_sl: int***    
    For some optimization methods, adjacent iterations have similar step lengths and descent directions. The optimal step length from the last iteration is an ideal test step length for the current iteration. *try_old_sl* indicates whether to use the optimal step length from the last iteration for step length estimation. 1 means yes; 0 means no.
  - ***max_m0: float***   
    The estimated maximum value of inversion parameters, which is useful for calculating the model update.
  - ***eps_scale: float***    
    The percentage of the model update relative to the model. If *eps_scale* > 0, the descent direction will be scaled to $eps\_scale*\frac{max\_m0}{max(p_k)}$; if *eps_scale* < 0, no scaling factor is applied to the descent direction. Instead, the user can apply the scaling factor to the gradient.
  - ***alpha_lb: float***  
    Lower limit of step length, suggested to be in the range 0 < $alpha\_lb$ < 1.
  - ***alpha_ub: float***  
    Upper limit of step length, suggested to be in the range 1 < $alpha\_ub$ < 10.
  - ***mod_limits: int***  
    Whether to check the bounds of model value. 1 means yes; 0 means no.
  - ***mod_lb: float***    
    Lower limit of model value, set according to the prior information of the model.
  - ***mod_ub: float***    
    Upper limit of model value, set according to the prior information of the model.

### 4.8 Inversion loop
  Assume that the inversion task has multiple stages, users need to specify the starting iteration, abort criteria, and so on.
  ```sh
  'Inversion loop'    :'comment',
  'abort'             :0, 
  'iter0'             :1,         # optim needs
  'niter_min'         :30,        # optim needs
  'niter_max'         :60,        # optim needs
  'pro'               :0.01,
  'pro2'              :0.01,
  'pro3'              :0.01,
  ```
  **Parameters**
  - ***abort: int***  
    The initial status for the inverse workflow. The value 0 is required when updating the configuration JSON file using the script **S0_prepare.py**; the value 1 indicates that one abort criterion has been met.
  - ***iter0: int***  
    The initial iteration number for the current inversion stage/phase.
  - ***niter_min: int***  
    The minimum number of iterations.
  - ***niter_max: int***  
    The maximum number of iterations, for future design.
  - ***pro: int***  
    The convergence threshold. The inversion will terminate when the relative reduction in the misfit function $\frac{E_{i-2}-E_{i}}{E_{i-2}}$ falls below this threshold, where $E_{i}$ means the misfit function at *i*th iteration.
  - ***pro2: int***  
    The slope convergence threshold. The inversion will terminate when the change rate in the misfit function slope $\frac{k_{i}}{k_{0}}$ falls below this threshold, where $k_{i}$ means the misfit function slope at *i*th iteration.
  - ***pro3: int***  
    The average convergence threshold. The inversion will terminate when the average value of the misfit function over nearby iterations $\frac{ave_{i}}{ave_{0}}$ falls below this threshold, where $ave_{i}$ means the average of misfit function around *i*th iteration.