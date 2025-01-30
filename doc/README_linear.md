# Guide to Linear Inverse Methods

>The *OPTIM* package provides a suite of common optimization algorithms designed to help users quickly build inversion workflows in the geophysical field. By leveraging this package, users can focus on the implementation of individual misfit functions and their first- and second-order derivatives. This optimization software comprises several popular descent direction modules and a line search strategy. This documentation will introduce two linear algorithms: Preconditioned Linear Conjugate Gradient (PLCG) and LSQR methods, for determining the descent direction. The optimization toolbox [optim](https://github.com/beyondTaoLei/optim) has been developed by Tao Lei (leit@sustech.edu.cn) and Wei Zhang (zhangwei@sustech.edu.cn). If you encounter any issues while using this software, please do not hesitate to contact us.
***

For convenience, the parameters of the linear inverse algorithms are organized into a configuration JSON file. This file includes the following key components: 1) program path module, 2) file path information, 3) descent direction module, 4) model update module, 5) inversion loop. Each of these components will be described in detail.

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
  We build the inverse workflow with a series of procedures, including the geophysical modules and the optimization modules, and each procedure connects the predecessors and successors with files. Therefore, it is necessary to specify the locations of the model, cost function, sensitivity matrix, and other relevant files. The contents of these files may be updated over iterations.
  ```sh
  'File path'         :'comment',
  'fmod'              :'model/inv/inv_tomo.vs',
  'fmodref'           :'model/inv/ref_tomo.vs',
  'fcost'             :'list/fcost.dat',
  'fdiff'             :'list/diff.dat',
  'fsensmat'          :'list/sensmat.npz'
  'fwgtm'             :'list/wgtm.npz',
  'fwgtd'             :'list/wgtd.npz',
  ```
  **Parameters**
  - ***fmod: input&output, string***  
    The path for storing the current model, where the data is saved in a 32-bit floating-point binary file, with a shape of (M, 1). The optimization module will read the model using *numpy.fromfile(fmod, np.float32)*. If the inversion target contains multiple parameters, these parameters need to be merged into a single vector for current package version.
  - ***fmodref: input, string***  
    Similar to ***fmod***, but the associated model is used as the reference model when the inversion seeks the smoothest model (***modeltype***=3).
  - ***fcost: output, string***  
    The path for storing the misfit function, where the data is saved in *ASCII* floating-point format.
  - ***fdiff: output, string***  
    The path for storing the measurement difference, where the data is saved in *ASCII* floating-point format, with a shape of (N, 1).
  - ***fsensmat: output, string***  
    The path for storing the sensitivity matrix, where the data is saved as a sparse matrix with a shape of (N, M).
  - ***fwgtm: optional input, string***  
    The path for storing the model weighting, where the data is saved as a compressed archive file in the format *.npz*. The data is a sparse diagonal matrix with a shape of (M, M).
  - ***fwgtd: optional input, string***  
    The path for storing the data weighting, where the data is saved as a compressed archive file in the format *.npz*. The data is a sparse diagonal matrix with a shape of (N, N).

### 3. Descent direction module
  For most optimization algorithms, users need to input few additional parameters to constrain the calculation of the descent direction. Here, we will introduce the basic parameters for linear optimization methods in the *OPTIM* software. You can choose either the Preconditioned Linear Conjugate Gradient (**PLCG**) or **LSQR** method to build your inverse workflow. In cases where the linear system is ill-posed or nearly singular, a regularization term can help stabilize the solutions by adding a controlled amount of bias, leading to more robust estimates. There are three types of regularization term to choose from, and more details can be found later.

#### 3.1 Preconditioned Linear Conjugate Gradient (PLCG) method
  The Linear Conjugate Gradient (LCG) method was proposed by Hestenes and Stiefel (1952), which is an iterative method mainly to solve linear equations, particularly those characterized by a positive definite matrix. This method allows us to approximately solve large linear systems at affordable iterations while the direct methods would take excessive computer time or storage, such as LU or Cholesky decomposition.
  ```sh
  'Descent part'      :'comment',
  'modeltype'         :3,                         
  'descent_type'      :1,                         
  'precond'           :1,                         
  'niter_inner_max'   :100,                      
  'lambda0'           :10000.0,                        
  ```
  **Parameters**
  - ***modeltype: int***  
    The desired model type one may look for, 1: the smoothest model, 2: the flattest model, 3: the smallest perturbation model.
  - ***descent_type: int***  
    The descent perturbation type when solving the linear equation. 1 indicates the absolute perturbation; 2 indicates the relative perturbation. 
  - ***precond: optional, int***  
    Whether to apply a preconditioning operator for the linear conjugate gradient method. 1 means yes; 0 means no.
  - ***niter_inner_max: int***  
    The maximum number of iterations for **PLCG** loop.
  - ***lambda0: float***  
    The regularization parameter is used to balance the data misfit and model misfit. The users typically need to make numerous attempts to identify the optimal parameters.

#### 3.2 LSQR method
  The Least-Squares with QR factorization (LSQR) represents another conjugate gradient algorithm designed for addressing large, ill-posed, and sparse least squares problems of any shape (Paige and Saunders, 1982). This algorithm is stable and efficient for solving the linear (or linearized) inversion problems, and can perform very few iterations to produce a satisfactory estimate. LSQR algorithm is similar in style to the linear conjugate gradient algorithm as applied to the linear inversion problems, so we implement these algorithms with the same interface in our proposed package.
  ```sh
  'Descent part'      :'comment',
  'modeltype'         :3,                         
  'descent_type'      :1,                         
  'niter_inner_max'   :100,                      
  'lambda0'           :10000.0,                        
  ```
  **Parameters**
  - ***modeltype: int***  
    The desired model type one may look for, 1: the smoothest model, 2: the flattest model, 3: the smallest perturbation model.
  - ***descent_type: int***  
    The descent perturbation type when solving the linear equation. 1 indicates the absolute perturbation; 2 indicates the relative perturbation. 
  - ***niter_inner_max: int***  
    The maximum number of iterations for **LSQR** loop.
  - ***lambda0: float***  
    The regularization parameter is used to balance the data misfit and model misfit. The users typically need to make numerous attempts to identify the optimal parameters.

### 4. model update module
  Once the descent direction $p_k$ is generated by the above linear optimization algorithms, the model can be updated either with a fixed step length or an optimal step length estimated by the parabolic fitting method. For linear inverse algorithms, we prefer to update the model with a fixed step length. Refer to the demo of first-arrival traveltime tomography (script *prob/update_model.py*) or surface wave inversion (script *prob/update_model.py*) in our package.
  
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