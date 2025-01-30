README for IFOS2D
=================

IFOS2D (Inversion of Full observed seismograms) 
is a 2D elastic full waveform inversion code. The inversion problem is solved 
by a conjugate gradient method. The gradients are computed in the time domain 
with the adjoint method. The forward modeling is also done in the time domain 
with a finite difference code on a standard staggered grid.


System requirements:
For a typical example of an inversion of shallow seismic surface waves with
- a model size of 400 grid points in lateral and 75 grid points in vertical
direction,
- 14000 time steps per forward modeling, 
- viscoelastic forward modeling with three relaxation mechanisms and
- 8 shots
we need
- at least 4 CPUs
- 720 MB memory
- around 4.6 MB storage capacities per iteration (in total around 800 MB for 173
iteration steps in this example)



Installation instructions can be found in the file INSTALL.


A detailed documentation is provided in manual_IFOS.pdf.


To get started with this code try to run the toy example that is described in
chapter 9 in the documentation (manual_IFOS.pdf).


If you use this code for your own research please cite at least one article
written by the developers of the package, e.g.
D. Köhn. Time domain 2D elastic full waveform tomography. PhD Thesis, Kiel
University, 2011.


Contact:
see file DEVELOPERS
