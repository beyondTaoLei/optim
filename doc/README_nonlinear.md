# Guide to Nonlinear Inverse Methods

>The *OPTIM* package provides a suite of common optimization algorithms designed to help users quickly build inversion workflows in the geophysical field. By leveraging this package, users can focus on the implementation of individual misfit functions and their first- and second-order derivatives. This optimization software comprises several popular descent direction modules and a line search strategy. This documentation will introduce four nonlinear algorithms: steepest descent (STD), nonlinear conjugate gradient (NLCG), l-BFGS, and truncated Newton (TRN) methods, for determining the descent direction. Additionally, the software employs a parabolic fitting method to estimate the update step length, a technique widely used in complex geophysical inversion workflows, such as full waveform inversion.
The optimization toolbox [optim](https://github.com/beyondTaoLei/optim) has been developed by Tao Lei (leit@sustech.edu.cn) and Wei Zhang (zhangwei@sustech.edu.cn). If you encounter any issues while using this software, please do not hesitate to contact us.
***

For convenience, the parameters of the nonlinear inverse algorithms are organized into a configuration JSON file. This file includes the following key components: 1) program path module, 2) file path information, 3) descent direction module, 4) model update module, 5) inversion loop. Each of these components will be described in detail.

### 1. Program path module
  Since we run all commands from the project root path, we need to define the location of the optimization toolbox. The following parameter is required:
  ```sh
  'Program path'      :'comment',
  'optimroot'         :'/path/to/software/optim',   # path to optim package
  ```
  **Parameters**
  - ***optimroot: string***   
  The root path of the optimization toolbox codes.

### 2. File path information
  Since we build the inverse workflow with a series of procedures, including the geophysical modules and the optimization modules, each procedure connects the predecessors and successors with files. Therefore, it is necessary to specify the locations of the model, cost function, gradient, etc. The contents of the following files may be updated over iterations.
  ```sh
  'File path'         :'comment',
  'fmod'              :'FD/mod.bin',
  'fcost'             :'FD/misfit.dat',
  'fgrad'             :'FD/grad.bin', 
  'fhess'             :'FD/hess.bin', 
  'fvctr'             :'FD/vctr.bin'
  ```
  **Parameters**
  - ***fmod: input&output, string***  
    The path for storing the current model, where the data is saved in a 32-bit floating-point binary file, with a shape of (M, 1). The optimization module will read the model using *numpy.fromfile(fmod, np.float32)*. If the inversion target contains multiple parameters, these parameters need to be merged into a single vector for current package version.
  - ***fcost: output, string***  
    The path for storing the misfit function, where the data is saved in *ASCII* floating-point format.
  - ***fgrad: output string***  
    The path for storing the current gradient, where the data is saved in a 32-bit floating-point binary file, with a shape of (M, 1).
  - ***fhess: optional output, string***  
    The path for storing the Hessian vector product *Hv*, where the data is saved in a 32-bit floating-point binary file, with a shape of (M, 1). This parameter is only necessary for the truncated Newton method.
  - ***fvctr: optional output, string***  
    The path for storing an arbitrary vector *v* in the truncated Newton method, where the data is saved in a 32-bit floating-point binary file, with a shape of (M, 1), This parameter is only necessary for the truncated Newton method.

### 3. Descent direction module
  For most optimization algorithms, users need to input few additional parameters to constrain the calculation of the descent direction. Here we will introduce the basic parameters for nonlinear optimization methods in the *OPTIM* software. You can choose any one optimization method to build your inverse workflow.

#### 3.1 Steepest descent direction (STD) method
  Since the descent direction is equal to the negative gradient, so no extra parameters are required.

#### 3.2 Nonlinear conjugate gradient (NLCG) method
  We prefer to use the formula proposed by Polak and Ribiere (1969) to construct the nonlinear conjugate gradient method in our study. Additionally, the Powell restart condition is added for fast convergence of inversion after generating the search descent direction (Powell, 1977).
  ```sh
  'NLCG'              :'comment',
  'powell'            :0                 
  ```
  **Parameters**
  - ***powell: int***  
    Whether to check the Powell restart condition, where 1 means the code will check this condition, while 0 means it will not.

  #### 3.3 l-BFGS method
  In general, the *l*-BFGS method utilizes the gradients and models of the most recent iterations to construct the descent direction, while discarding information from earlier iterations to save storage.
  ```sh
  'l-BFGS'            :'comment',
  'l'                 :5
  ```
  **Parameters**
  - ***l: int***  
    The number of the most recent correction pairs $\{s_k, y_k\}$. Users can specify the number *l* according to the degree of complexity of the inverse problem. It is suggested to be less than 10.
  
  #### 3.4 Truncated Newton (TRN) method
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
  'eta'               :0.9,               
  'eta_update'        :0                     
  ```
  **Parameters**
  - ***conv_CG: int***  
    The initial status for the inner loop. The value 0 is required when you update the configuration JSON file using the script **S0_prepare.py**.
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

### 4. Model update module
  Once the descent direction $p_k$ is generated by the above nonlinear optimization algorithms, line search strategy can help to find the optimal step length $\alpha_k$ for the model updating. We choose the parabolic fitting method to achieve this goal, and the relevant parameters are as follows:
  ```sh
  'Steplength'        :'comment',
  'try_old_sl'        :1,                     
  'max_m0'            :0.25,                
  'eps_scale'         :0.1,                 
  'alpha_lb'          :0.01,                   
  'alpha_ub'          :5,                    
  'mod_limits'        :0,                     
  'mod_lb'            :1000,                  
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

### 5. Inversion loop
  Assume that the inversion task has multiple stages, and users need to specify the starting iteration for the current stage.
  ```sh
  'Inversion loop'    :'comment',
  'iter0'             :1,                     
  'niter_max'         :30
  ```
  **Parameters**
  - ***iter0: int***  
    The initial iteration number at the current inversion stage/phase.
  - ***niter_max: int***  
    The maximum number of iterations, for future design.