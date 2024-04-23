# CONTENTS

The OPTIM package provides some common optimization algorithms, aiming to help the potential users to quickly build the inversion workflow in geophysical field, while the users can focus on implementing the individual misfit, and its first- and second-order derivatives. This optimization software consists of several popular descent direction modules, and line search strategy. there are four nonlinear algorithms: steepest descent (STD), nonlinear conjugate gradient (NLCG), l-BFGS and truncated Newton (TRN) methods, and two linear algorithms: linear conjugate gradient (LCG) and LSQR methods to generate the descent direction in this package. It also provide the option of incorporating the regularization terms to constrain the model construction. Furthermore, the software employs the parabolic fitting method to estimate the update step length, a technique widely utilized in complex geophysical inversion workflows like full waveform inversion.

This software works on few common open-source Python libraries, including Numpy and Scipy, which facilitate software installation for the potential users. To enhance user experience, we offer multiple demonstrations covering various inverse problems, along with comprehensive technical documentation. The overall software structure is illustrated in Figure \ref{fig:optim}.

&#x20;To successfully execute the scripts in the demonstrations, please install the additional Python libraries as follows:

*   **First arrival traveltime tomography**: pandas, matplotlib

*   **Full waveform inversion (FWI)**: pandas, matplotlib, mpi4py, obspy

*   **Surface wave inversion**: pandas, matplotlib, mpi4py, disba, and DisbaTomo&#x20;

# PRINCIPLES

Considering the runtime of optimization module is much smaller than that of geophysical module, we choose interpreted language Python to implement our optimization module. Python is known for being portable, productive, open source, and easy to read. What's more, there are large libraries available for users, such as Numpy, Scipy and Matplotlib, etc. Object-oriented programming has been popular in recent decades, and its inheritance and polymorphism prove useful in most cases, but they can also add unnecessary complexity to the code structure. From the creator's point of view, the above features are conducive to the modular management of software. However, the users often need to thoroughly track the inheritance relationship of parent and child classes in multiple files to understand the data flow in memory. In this paper, the design of proposed software OPTIM follows two principles:

*   **Writing code as the underlying mathematics.**
    The implementation of each algorithm is mostly completed in one file, avoiding the logical complexity caused by function jumping as much as possible, and it also helps the users understand and further improve the corresponding algorithm. As Wikipedia "Unix philosophy" emphasizes, the proposed codes are expected to be easily maintained and repurposed by future developers rather than its creators.

*   **Assemble the inverse workflow by writing short programs.**
    All programs start to work from a shared configuration which saves the necessary parameters and file paths. The output files from the former steps will be entered into the latter step. "Do one thing and do it well" will keep the whole workflow simple and clear.

