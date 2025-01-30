# Surface wave inversion Documentation

>Surface wave inversion is a geophysical method used to infer subsurface properties by analyzing the dispersion characteristics of surface waves, particularly Rayleigh and Love waves. This technique is widely applied in engineering, environmental studies, and earthquake seismology to obtain shear wave velocity profiles. In this approach, the primary observational data are the phase (or group) velocities of surface waves, which vary with frequency due to the layered structure of the Earth's subsurface. By measuring the phase velocity at different frequencies, one can obtain dispersion curves that are indicative of the subsurface velocity structure. The inversion process involves fitting observed dispersion curves with theoretical ones to derive a velocity model that best represents the subsurface layers. This model provides crucial information about the elastic properties and layering of the subsurface materials.
This documentation aims to help users understand our implementation process for solving **medium-scale** geophysical inverse problems. For large-scale inverse problems, refer to the full waveform inversion demo (fwi) included in our packages, which will guide you on how to organize the optimization workflow and submit jobs to an HPC cluster. The optimization toolbox [optim](https://github.com/beyondTaoLei/optim) has been developed by Tao Lei (leit@sustech.edu.cn) and Wei Zhang (zhangwei@sustech.edu.cn). If you encounter any issues while using this software, please do not hesitate to contact us.
***

## 0. Demo structure
The root directory contains several key components essential for understanding and utilizing the provided resources. It includes this documentation, the workflow example (**case1**), the collection of main functions (**prob**), and some drawing scripts located in **plot**. 

Most parameters will be set using the script **S0_prepare.py** for linear inverse algorithms or **S0_prepareNL.py** for nonlinear inverse algorithms within the workflow example(**case1**). The following directory tree presents the main structure for solving first arrival traveltime tomography problem:

```
|-- README.md                   >>> this documentation
|-- case1                       >>> workflow example
|   |-- P0_prepare_files.py
|   |-- P1_generate_data.py
|   |-- S0_prepare.py           >>> JSON for linear inverse algo.
|   |-- S1_optim_LSQR.py
|   |-- S1_optim_PLCG.py
|   |-- S0_prepareNL.py         >>> JSON for nonlinear inverse algo.
|   |-- S1_optim_STD.py
|   |-- S1_optim_NLCG.py
|   |-- S1_optim_LBFGS.py
|   |-- S1_optim_TRN.py
|   `-- S2_run.py
|-- plot                        >>> figure presentation
|   |-- plot_mod.py
|   |-- plot_disp.py
|   |-- plot_misfit.py
`-- prob                        >>> main functions used in workflow
    |-- forward.py
    |-- kernel.py
    |-- calc_diff_sensmat.py
    |-- media_fd2inv.py
    |-- my_misfit.py
    |-- my_grad.py
    |-- my_Hv.py
    |-- update_model.py
    |-- update_optim.py
    |-- check_stopping.py
    `-- util.py

      Figure 1. The structure of surface wave inversion demo
```

## 1. Configure environment

Our surface wave inversion utilizes the computational capabilities of the softwares *disba* and *DisbaTomo* for calculating the dispersion and phase velocity sensitivity kernels to shear wave velocity. Follow the steps below to set up the required environment for the successful installation and compilation of our software.

- **Prerequisites**
  Please ensure you have the following prerequisites in place:
*Python* (version >= 3.8)
*pip* (Python package installer)
*sudo* privileges may be required for certain installations.

- **Install DisbaTomo**
To install DisbaTomo, refer to the installation commands provided in the [DisbaTomo GitHub repository](https://github.com/pan3rock/DisbaTomo) under the sections **With sudo permission** and **Compilation**. Follow the detailed instructions to clone the repository and compile the library. After completing the installation and compilation, verify that *DisbaTomo* is correctly installed by running a simple test script provided in their documentation.

- **Copy the compiled library**: copy the compiled shared object file, as *dpcxx.cpython-[python version]-[processor architecture]-linux-gnu.so*, to the *surfacewave* demo directory of our software. for example:
  ```sh
  cp /path/to/your/DisbaTomo/software/python/dpcxx.cpython-312-x86_64-linux-gnu.so \
     /path/to/your/optim/software/demo/surfacewave/prob/
  ```

By following these steps, you will have set up the necessary environment to utilize our software effectively.

Note: For detailed installation instructions and troubleshooting, refer to the official documentation of *DisbaTomo*.

## 2. Prepare observed data

- Specify the path to the tomography demo and name your project in the terminal
  ```sh
  demoroot=/path/to/software/optim/demo/surfacewave
  mkdir prj
  cd prj
  ```
  Here, **demoroot** indicates the root path of the example codes. Specifying this path in the terminal will facilitate copying and using the subsequent codes for study.

- Modify the program path variables in the script **${demoroot}/case1/S0_prepare.py**  
  ```sh
  'Program path'      :'comment',
  'optimroot'         :'/path/to/software/optim',
  ```

- Write the global parameters for the whole project
  ```sh
  python ${demoroot}/case1/S0_prepare.py optim
  ```
  :bell: Output: **optim/optim.json**
  We place all optimization parameters in the script **S0_prepare.py**. Each program in the inverse workflow will extract the necessary parameters from the output JSON file. The input parameter **optim** in the command specifies the output folder to store the global parameters.

- Copy the models from the *optim* package into the current project directory, distribute the true model to a specific location, and display the model as
  ```sh
  cp -r ${demoroot}/case1/model .
  cp model/true/true.vp model/inv/inv.vp
  cp model/true/true.vs model/inv/inv.vs
  cp model/true/true.rho model/inv/inv.rho
  python ${demoroot}/plot/plot_mod.py optim model/inv/inv.vs
  ```
  | Target          | :bell: Output |
  | ------          | ------ |
  | vp model        | **model/inv/inv.vp** |
  | vs model        | **model/inv/inv.vs** |
  | density model   | **model/inv/inv.rho** |

- Prepare the geometry for the current study, which includes 129 stations located on the surface, and the measurement periods range from 3 to 59 s with an interval of 2 s. 
  ```sh
  python ${demoroot}/case1/P0_prepare_files.py optim
  ``` 
  | Target          | :bell: Output |
  | ------          | ------ |
  | station         | **list/stn.csv** |
  | period          | **list/periods.csv** |

- Generate the observed data for all stations
  ```sh
  # first generate the commands
  python ${demoroot}/case1/P1_generate_data.py optim
  # and then execute with MPI processors
  python ${demoroot}/case1/S2_run.py optim
  ```
  | Target              | :bell: Output |
  | ------              | ------ |
  | job list            | **job/job.csv** |
  | observed dispersion | **obs/vmap_M0.bin** |
  
  where the first script will generate the forward modeling commands and store them in the file **job/job.csv**, while the second script will execute these commands and generate the observed dispersion. The number after the letter *M* represents the mode of the dispersion data, 0 means the fundamental mode, 1 means the first higher mode and so on.

- Display the observed dispersion data
  ```sh
  python ${demoroot}/plot/plot_disp.py optim obs/vmap_M0.bin
  ```

## 3. Run linear inverse workflow

- Prepare the initial model to start the inversion process
  ```sh
  cp model/init/init.vp model/inv/inv.vp
  cp model/init/init.vs model/inv/inv.vs
  cp model/init/init.rho model/inv/inv.rho
  ```
  | Target          | :bell: Output |
  | ------          | ------ |
  | vp model        | **model/inv/inv.vp** |
  | vs model        | **model/inv/inv.vs** |
  | density model   | **model/inv/inv.rho** |
  
- Generate the LSQR workflow from iteration 1 to 20, and then execute the workflow
  ```sh
  # prepare the global parameters
  python ${demoroot}/case1/S0_prepare.py LSQR
  # generate the commands
  python ${demoroot}/case1/S1_optim_LSQR.py LSQR 1 20
  # execute with MPI processors
  python ${demoroot}/case1/S2_run.py LSQR
  ```
  | Target              | :bell: Output |
  | ------              | ------ |
  | job list            | **job/job.csv** |
  | evaluated model     | **LSQR/iter#/modg.bin** |
  | data misfit         | **LSQR/iter#/fcost.dat** |

  where the file **job.csv** contains the commands for the inverse workflow, **modg.bin** represents the evaluated model at the **#**-th iteration, and **fcost.dat** represents the total data misfit.

- Display the inverse results using the LSQR approach
  ```sh
  # plot model
  python ${demoroot}/plot/plot_mod.py LSQR LSQR/iter20/modg.bin
  # plot data
  cp syn/vmap_M0.bin syn/vmap_M0_LSQR20.bin
  python ${demoroot}/plot/plot_disp.py LSQR syn/vmap_M0_LSQR20.bin
  # plot misfit func.
  python ${demoroot}/plot/plot_misfit.py LSQR
  ```

- Test PLCG optimization method in the same way (_optional_)
  ```sh
  # initialize the workflow
  cp model/init/init.vp model/inv/inv.vp
  cp model/init/init.vs model/inv/inv.vs
  cp model/init/init.rho model/inv/inv.rho
  # prepare the global parameters
  python ${demoroot}/case1/S0_prepare.py PLCG
  # generate the commands
  python ${demoroot}/case1/S1_optim_PLCG.py PLCG 1 20
  # execute with MPI processors
  python ${demoroot}/case1/S2_run.py PLCG
  # plot model
  python ${demoroot}/plot/plot_mod.py PLCG PLCG/iter20/modg.bin
  # plot data
  cp syn/vmap_M0.bin syn/vmap_M0_PLCG20.bin
  python ${demoroot}/plot/plot_disp.py PLCG syn/vmap_M0_PLCG20.bin
  # plot misfit func.
  python ${demoroot}/plot/plot_misfit.py PLCG
  ```
  where the parameter **PLCG** following the commands specifies a folder used to store the results from each iteration.

## 4. Set json Parameters for linear inverse algorithms

In this section, we will discribe how to set the parameters in the configuration JSON file. The final file should include the following key components: 1) program path module, 2) file path information, 3) forward modeling module, 4) descent direction module, 5) model update module, 6) inversion loop. We will go over each of these components in detail.

### 4.1 Program path module
  Since we run all commands from the project root path, we need to define the location of the optimization toolbox. The following parameter is required:
  ```sh
  'Program path'      :'comment',
  'optimroot'         :'/path/to/software/optim',
  ```
  **Parameters**
  - ***optimroot: string*** 
  The root path to the optimization toolbox codes. 

### 4.2 File path information
  We build the inverse workflow with a series of procedures, including the geophysical modules and the optimization modules, and each procedure connects the predecessors and successors with files. Therefore, it is necessary to specify the locations of the model, cost function, sensitivity matrix, and other relevant files. The contents of these files may be updated over iterations.
  ```sh
  'File path'         :'comment',
  'fstn'              :'list/stn.csv',
  'fperiods'          :'list/periods.csv',
  'fmodfd'            :'model/inv/inv.vs',
  'fmodfdref'         :'model/inv/ref.vs',
  'fmod'              :'model/inv/inv_tomo.vs',   # optim needs
  'fmodref'           :'model/inv/ref_tomo.vs',   # optim needs
  'fcost'             :'list/fcost.dat',          # optim needs
  'fdiff'             :'list/diff.dat',           # optim needs
  'fsensmat'          :'list/sensmat.npz'         # optim needs
  'fwgtm'             :'list/wgtm.npz',           # optim needs
  'fwgtd'             :'list/wgtd.npz',           # optim needs
  ```
  **Parameters**
  - ***fstn: input, string***
    The path for storing the station information, where the data is saved in *csv* format. The key 'station' in the header indicates the station name, while 'x', 'y' and 'z' indicate the coordinates of the station, and 'ix', 'iy' and 'iz' indicate the grid index of the station.
  - ***fperiods: input, string***
    The path for storing the measurement periods, where the data is saved in *csv* format. The key 'phase' in the header indicates the period index, while 'T0' indicates the measurement period. Note that the period must be increasing.
  - ***fmodfd: input&output, string***
    The path for storing the current vs model, where the data is saved in a 32-bit floating-point binary file. The input data is arranged in *column-major* order, with the coordinate origin at the top-left corner. The x-axis extends to the right, and the y-axis extends downward. Note that, the associated vp and density models have also to be provided, named with the suffixes .vp and .rho, respectively.
  - ***fmodfdref: optional input, string***
    Similar to ***fmodfd***, but the associated model is used as the reference model when the inversion seeks the smoothest model (***modeltype***=3).
  - ***fmod: output, string***
    The path for storing the resampled output of the model is within the file ***fmodfd***, with a shape of (M, 1).
  - ***fmodref: output, string***
    The path for storing the resampled output of the model is within the file ***fmodfdref***, with a shape of (M, 1).
  - ***fcost: output, string***
    The path for storing the misfit function, where the data is saved in *ASCII* floating-point type.
  - ***fdiff: output, string***
    The path for storing the measurement difference, where the data is saved in *ASCII* floating-point format, with a shape of (N, 1).
  - ***fsensmat: output, string***
    The path for storing the sensitivity matrix, where the data is saved as a sparse matrix with a shape of (N, M).
  - ***fwgtm: optional input, string***
    The path for storing the model weighting, where the data is saved as a compressed archive file in the format *.npz*. The data is a sparse diagonal matrix with a shape of (M, M).
  - ***fwgtd: optional input, string***
    The path for storing the data weighting, where the data is saved as a compressed archive file in the format *.npz*. The data is a sparse diagonal matrix with a shape of (N, N).

### 4.3 Forward modeling module
  Next, we will detail the parameters specifically for travel time and sensitivity calculations. 
  ```sh
  'Forward part'      :'comment',
  'mpiexep'           :'mpirun',
  'algo'              :'dunkin',
  'wavetype'          :'rayleigh',
  'maxmode'           :0,
  'ncores'            :8,
  'nxyzfd'            :[1, 1300, 400],
  'oxyzfd'            :[0.0, 0.0, 0.0],
  'dxyzfd'            :[500.0, 500.0, 500.0],
  'incxyz'            :[1, 10, 1],
  ```
  **Parameters**
  - ***mpiexep: int***
    The path to Python MPI parallel launcher.
  - ***algo: int***
    The algorithm to use for computation of Rayleigh-wave dispersion, *dunkin* and *fast-delta* is supported by the package *disba*.
  - ***wavetype: int***
    The wave type, one can specify the wave type with *love* or *rayleigh*.
  - ***maxmode: int***
    The maximum modes of the dispersion measurements, 0 means the fundamental mode, 1 means the first higher mode and so on.
  - ***ncores: int***
    The number of cores for the calculation of dispersion curve and sensitivity kernel, and we implement the calculation with MPI processors.
  - ***nxyzfd: list of int***
    The number of grid points in the x-, y- and z-directions for forward modeling, respectively.
  - ***oxyzfd: list of float***
    The original coordinates of the input model in the x-, y- and z-directions for forward modeling, respectively.
  - ***dxyzfd: list of float***
    The grid spacing of the input model in the x-, y- and z-directions for forward modeling, respectively.
  - ***incxyz: list of int***
    The sampling interval number in the x-, y- and z-directions, respectively, for saving the sensitivity matrix.

### 4.4 Descent direction module
  For most optimization algorithms, users need to input few additional parameters to constrain the calculation of the descent direction. Here, we will introduce the basic parameters for linear optimization methods in the *OPTIM* software. You can select the **LSQR** or Preconditioned Linear Conjugate Gradient (**PLCG**) method to build your inverse workflow. In cases where the linear system is ill-posed or nearly singular, a regularization term can help stabilize the solutions by adding a controlled amount of bias, leading to more robust estimates. There are three types of regularization term to choose from, and more details can be found later.
  ```sh
  'Descent part'      :'comment',
  'modeltype'         :3,             # optim needs                         
  'descent_type'      :1,             # optim needs            
  'precond'           :0,             # optim needs            
  'niter_inner_max'   :100,           # optim needs           
  'lambda0'           :6.0e-3,        # optim needs                 
  ```
  **Parameters**
  - ***modeltype: int***
    The desired model type one may look for, 1: the smoothest model, 2: the flattest model, 3: the smallest perturbation model.
  - ***descent_type: int***
    The descent perturbation type when solving the linear equation. 1 indicates the absolute perturbation; 2 indicates the relative perturbation. 
  - ***precond: optional, int***
    Whether to apply a preconditioning operator for the linear conjugate gradient method. 1 means yes; 0 means no.
  - ***niter_inner_max: int***
    The maximum number of iterations for **PLCG** or **LSQR** loop.
  - ***lambda0: float***
    The regularization parameter is used to balance the data misfit and model misfit. The users typically need to make numerous attempts to identify the optimal parameters.

### 4.5 Model update module
  Once the descent direction $p_k$ is generated by the above linear optimization algorithms, we will update the model by a fixed percentage. The relevant parameters are as follows,
  ```sh
  'Update part'       :'comment',
  'eps_scale'         :0.03,                      
  'smooth_size'       :[0.0, 30000.0, 9000.0],            
  'mod_limits'        :1,                         
  'mod_lb'            :2500,                      
  'mod_ub'            :5000,                      
  ```
  **Parameters**
  - ***eps_scale: float***
    The percentage by which the model is updated relative to its current state.
  - ***smooth_size: list of float*** 
    The smoothing size in the x, y and z directions for the descent direction, measured in meters.
  - ***mod_limits: int***
    Whether to check the bounds of model value. 1 means yes; 0 means no.
  - ***mod_lb: float***  
    Lower limit of model value, set according to the prior information of the model.
  - ***mod_ub: float***  
    Upper limit of model value, set according to the prior information of the model.

### 4.6 Inversion loop
  Assume that the inversion task has multiple stages, users need to specify the starting iteration, abort criteria, and so on.
  ```sh
  'Inversion loop'    :'comment',
  'abort'             :0, 
  'iter0'             :1,         # optim needs
  'niter_min'         :20,        # optim needs
  'niter_max'         :30,        # optim needs
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

## 5. Run nonlinear inverse workflow

- Prepare the initial model to start the inversion process
  ```sh
  cp model/init/init.vp model/inv/inv.vp
  cp model/init/init.vs model/inv/inv.vs
  cp model/init/init.rho model/inv/inv.rho
  ```
  | Target          | :bell: Output |
  | ------          | ------ |
  | vp model        | **model/inv/inv.vp** |
  | vs model        | **model/inv/inv.vs** |
  | density model   | **model/inv/inv.rho** |
  
- Generate the LBFGS workflow from iteration 1 to 30, and then execute the workflow
  ```sh
  # prepare the global parameters
  python ${demoroot}/case1/S0_prepareNL.py LBFGS
  # generate the commands
  python ${demoroot}/case1/S1_optim_LBFGS.py LBFGS 1 30
  # execute with MPI processors
  python ${demoroot}/case1/S2_run.py LBFGS
  ```
  | Target              | :bell: Output |
  | ------              | ------ |
  | job list            | **job/job.csv** |
  | evaluated model     | **LBFGS/iter#/mod.bin** |
  | data misfit         | **LBFGS/iter#/fcost.dat** |
  | gradient            | **LBFGS/iter#/grad.bin** |

  where the file **job.csv** contains the commands for the inverse workflow, **mod.bin** represents the evaluated model at the **#**-th iteration, and **fcost.dat** and **grad.bin** represents the total data misfit and the model gradient.

- Display the inverse results using the LBFGS approach
  ```sh
  # plot model
  python ${demoroot}/plot/plot_mod.py LBFGS LBFGS/iter30/mod.bin
  # plot data
  cp syn/vmap_M0.bin syn/vmap_M0_LBFGS30.bin
  python ${demoroot}/plot/plot_disp.py LBFGS syn/vmap_M0_LBFGS30.bin
  # plot misfit func.
  python ${demoroot}/plot/plot_misfit.py LBFGS
  ```

To avoid repetition, we will not detail the steps for nonlinear inversion again. To achieve great inversion results, The recommended number of iterations for various inversions are as follows:
  | Optimization    |workflow script    |Iteration No.  |
  | ------          | ------            |  ------       |
  |LSQR             |S1_optim_LSQR.py   | *20*          |
  |PLCG             |S1_optim_PLCG.py   | *20*          |
  |STD              |S1_optim_STD.py    | *50*          |
  |NLCG             |S1_optim_NLCG.py   | *50*          |
  |LBFGS            |S1_optim_LBFGS.py  | *30*          |
  |TRN              |S1_optim_TRN.py    | *20*          |

We recommend that users refer to the LBFGS inversion process when testing the aforementioned nonlinear inversion algorithms.

## 6. Set json Parameters for nonlinear inverse algorithms

In this section, we will discribe how to set the parameters in the configuration JSON file. The final file should include the following key components: 1) program path module, 2) file path information, 3) forward modeling module, 4) gradient calculation module, 5) descent direction module, 6) model update module, 7) inversion loop. We will go over each of these components in detail.

### 6.1 Program path module
  Since we run all commands from the project root path, we need to define the location of the optimization toolbox. The following parameter is required:
  ```sh
  'Program path'      :'comment',
  'optimroot'         :'/path/to/software/optim',
  ```
  **Parameters**
  - ***optimroot: string*** 
  The root path to the optimization toolbox codes. 

### 6.2 File path information
  We build the inverse workflow with a series of procedures, including the geophysical modules and the optimization modules, and each procedure connects the predecessors and successors with files. Therefore, it is necessary to specify the locations of the model, cost function, sensitivity matrix, and other relevant files. The contents of these files may be updated over iterations.
  ```sh
  'File path'         :'comment',
  'fstn'              :'list/stn.csv',
  'fperiods'          :'list/periods.csv',
  'fmodfd'            :'model/inv/inv.vs',
  'fmod'              :'model/inv/inv.vs',    # optim needs
  'fcost'             :'list/fcost.dat',      # optim needs
  'fgrad'             :'list/grad.bin',       # optim needs
  'fhess'             :'list/hess.bin',       # optim needs
  'fvctr'             :'list/vctr.bin',       # optim needs
  ```
  **Parameters**
  - ***fstn: input, string***
    The path for storing the station information, where the data is saved in *csv* format. The key 'station' in the header indicates the station name, while 'x', 'y' and 'z' indicate the coordinates of the station, and 'ix', 'iy' and 'iz' indicate the grid index of the station.
  - ***fperiods: input, string***
    The path for storing the measurement periods, where the data is saved in *csv* format.
  - ***fmodfd: input&output, string***
    The path for storing the current vs model, where the data is saved in a 32-bit floating-point binary file. The input data is arranged in *column-major* order, with the coordinate origin at the top-left corner. The x-axis extends to the right, and the y-axis extends downward. Note that, the associated vp and density models have also to be provided, named with the suffixes .vp and .rho, respectively.
  - ***fmod: output, string***
    The path for storing the resampled output of the model is within the file ***fmodfd***, with a shape of (M, 1). For this study, we do not resample the grids during calculating the gradient.
  - ***fcost: output, string***
    The path for storing the misfit function, where the data is saved in *ASCII* floating-point type.
  - ***fgrad: output, string***
    The path for storing the gradient, where the data is saved in a 32-bit floating-point binary file.
  - ***fhess: optional output, string***
    The path for storing the Hessian vector product, where the data is saved in a 32-bit floating-point binary file, and is only necessary for the truncated Newton method.
  - ***fvctr: optional output, string***
    The path for storing an arbitrary vector in the LCG algorithm, where the data is saved in a 32-bit floating-point binary file, and is only necessary for the truncated Newton method.

### 6.3 Forward modeling module
  Next, we will detail the parameters specifically for travel time and sensitivity calculations. 
  ```sh
  'Forward part'      :'comment',
  'mpiexep'           :'mpirun',
  'algo'              :'dunkin',
  'wavetype'          :'rayleigh',
  'maxmode'           :0,
  'ncores'            :8,
  'nxyzfd'            :[1, 1300, 400],
  'oxyzfd'            :[0.0, 0.0, 0.0],
  'dxyzfd'            :[500.0, 500.0, 500.0],
  ```
  **Parameters**
  - ***mpiexep: int***
    The path to Python MPI parallel launcher.
  - ***algo: int***
    The algorithm to use for computation of Rayleigh-wave dispersion, *dunkin* and *fast-delta* is supported by the package *disba*.
  - ***wavetype: int***
    The wave type, one can specify the wave type with *love* or *rayleigh*.
  - ***maxmode: int***
    The maximum modes of the dispersion measurements, 0 means the fundamental mode, 1 means the first higher mode and so on.
  - ***ncores: int***
    The number of cores for the calculation of dispersion curve and sensitivity kernel, and we implement the calculation with MPI processors.
  - ***nxyzfd: list of int***
    The number of grid points in the x-, y- and z-directions for forward modeling, respectively.
  - ***oxyzfd: list of float***
    The original coordinates of the input model in the x-, y- and z-directions for forward modeling, respectively.
  - ***dxyzfd: list of float***
    The grid spacing of the input model in the x-, y- and z-directions for forward modeling, respectively.

### 6.4 Gradient calculation module
  Next, we will detail the parameters specifically for gradient calculations, and the associated parameters are as
  ```sh
  'Gradient part'     :'comment',
  'smooth_size'       :[0.0, 30000.0, 9000.0],    
  'precond'           :1,                         
  'EPSILON'           :0.05,
  ```
  **Parameters**
  - ***smooth_size: list of float*** 
    The smoothing size in the x, y and z directions for the descent direction, measured in meters.
  - ***precond: int*** 
    Whether to apply a preconditioning operator for the gradient. 1 means yes; 0 means no.
  - ***EPSILON: float*** 
    The water level to stabilize the preconditioning operator, suggested in 0~0.05.

### 6.5 Descent direction module
  For most optimization algorithms, users need to input few additional parameters to constrain the calculation of the descent direction. Here we will introduce the basic parameters for nonlinear optimization methods in the *OPTIM* software. You can choose any one optimization method to build your inverse workflow.

#### 6.5.1 Steepest descent direction (STD) method
  Since the descent direction is equal to the negative gradient, so no extra parameters are required.

#### 6.5.2 Nonlinear conjugate gradient (NLCG) method
  We prefer to use the formula proposed by Polak and Ribiere (1969) to construct the nonlinear conjugate gradient method in our study. Additionally, the Powell restart condition is added for fast convergence of inversion after generating the search descent direction (Powell, 1977).
  ```sh
  'NLCG'              :'comment',
  'powell'            :0                 
  ```
  **Parameters**
  - ***powell: int***
    Whether to check the Powell restart condition, where 1 means the code will check this condition, while 0 means it will not.

  #### 6.5.3 l-BFGS method
  In general, the *l*-BFGS method utilizes the gradients and models of the most recent iterations to construct the descent direction, while discarding information from earlier iterations to save storage.
  ```sh
  'l-BFGS'            :'comment',
  'l'                 :5
  ```
  **Parameters**
  - ***l: int***
    The number of the most recent correction pairs $\{s_k, y_k\}$. Users can specify the number *l* according to the degree of complexity of the inverse problem. It is suggested to be less than 10.
  
  #### 6.5.4 Truncated Newton (TRN) method
  When performing the truncated Newton method to update the model, two loops are involved. The outer loop aims to find the local minimum of the misfit function, while the inner loop aims to calculate the Newton direction by solving the Newton equation. We apply the linear conjugate gradient (LCG) method to generate the descent direction.

  To ensure that the solution is a descent direction and that the fast convergence rate of the pure Newton method is preserved, we add several tests when generating the Newton direction. These include singularity test, (strong) negative curvature test, and truncated test (Xie and Schlick, 1999). If any of these test conditions are met, we terminate the LCG algorithm. There are 11 parameters to constrain the Newton descent direction. 
  ```sh
  'TRN'               :'comment',
  'conv_CG'           :0,                 
  'niter_max_CG'      :12,                
  'flag_singularity'  :1,                 
  'flag_negative'     :2,                 
  'flag_truncation'   :1,                 
  'flag_min_res'      :1,                 
  'pro_sing'          :1.0E-10,           
  'pro_curv'          :1.0E-10,           
  'pro_trun'          :0.5,               
  'eta_update'        :0                     
  'eta'               :0.0,               
  ```
  **Parameters**
  - ***conv_CG: int***
    The initial status for the inner loop. The value 0 is required when you update the configuration JSON file using the script **S0_prepareNL.py**.
  - ***niter_max_CG: int***
    The maximum number of inner LCG loop, suggested to be about 10-12 for large-scale nonlinear inverse problems.
  - ***flag_singularity: int*** 
    Whether to perform the singularity test. Set to 1 for yes, and 0 for no.
  - ***flag_negative: int*** 
    Whether to perform the negative curvature test. The values can be 1 to perform the negative curvature test, 2 to perform the strong negative curvature test, or 0 to not perform any curvature test.
  - ***flag_truncation: int*** 
    Whether to perform the truncated test. The values can be 1 to perform the truncation test based on residuals, or 0 to not perform the test.
  - ***flag_min_res: int*** 
    Whether to replace the descent direction by the one when the residual of the Newton equation is at a minimum. The values can be 1 to perform the replacement, or 0 to not perform the replacement.
  - ***pro_sing: float*** 
    The threshold $\varsigma$ for the singularity test. The suggested value is 1.0E-10.
  - ***pro_curv: float***
    The threshold $\delta$ for the negative curvature test (*flag_negative* = 1). The suggested value is 1.0E-10.
  - ***pro_trun: float*** 
    The threshold $c_q$ for the truncation test. The suggested value is 0.5.
  - ***eta_update: int*** 
   Whether to perform and update the forcing term by Eisenstat and Walker (1994) before the inner LCG loop. The values can be 1 to check and update the *eta* during the inversion cycle, or 0 to not check this forcing term.
  - ***eta: float***
    Useful when *eta_update* is set to 1. The initial value is suggested to be 0.0~1.0.

### 6.6 Model update module
  Once the descent direction $p_k$ is generated by the above nonlinear optimization algorithms, line search strategy can help to find the optimal step length $\alpha_k$ for the model updating. We choose the parabolic fitting method to achieve this goal, and the relevant parameters are as follows:
  ```sh
  'Steplength'        :'comment',
  'try_old_sl'        :1,                     
  'max_m0'            :4700.0,                
  'eps_scale'         :0.03,                 
  'alpha_lb'          :0.1,                   
  'alpha_ub'          :10,                    
  'mod_limits'        :1,                     
  'mod_lb'            :2500,                  
  'mod_ub'            :5000                  
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

### 6.7 Inversion loop
  Assume that the inversion task has multiple stages, users need to specify the starting iteration, abort criteria, and so on.
  ```sh
  'Inversion loop'    :'comment',
  'abort'             :0, 
  'iter0'             :1,         # optim needs
  'niter_min'         :30,        # optim needs
  'niter_max'         :40,        # optim needs
  'pro'               :0.01,
  'pro2'              :0.01,
  'pro3'              :0.01,
  ```
  **Parameters**
  - ***abort: int***
    The initial status for the inverse workflow. The value 0 is required when updating the configuration JSON file using the script **S0_prepareNL.py**; the value 1 indicates that one abort criterion has been met.
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