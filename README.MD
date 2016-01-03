********************************************************************

	Diffusion Low-Field MRI Monte Carlo (DMMC)
	 Simulations on Multi-GPU Systems
   Author: Trevor Vincent (tvincent@cita.utoronto.ca)
          Last Updated: September 17, 2013

********************************************************************


Intro
-----

The general ideas behind the DMMC simulations are explained in my Honour's Thesis (email me for it).
Please look there for explanations of the algorithms and deeper analysis of the underlying backend (especially look at the appendix). 
This Multi-GPU code supports general regions with multiple compartments
having different D's, T2's, T1's and P's and any general magnetic field sequences (i.e. containing oscillating RF fields,
concomitant terms etc). Since speed is very important in the GPU version, the code is split into a general version which
makes no assumptions and solves the Bloch equations and a version which makes the high-field approximation and just 
updates the phase. My DMMC code has a frontend which is decently easy to use and allows for very complicated geometries and sequences
to be implemented in very few lines of code (the difference between a FCC lattice of spheres and a single sphere is just 13 lines of
frontend code!).

In this readme I will go over how to compile and run theGPU version. Then I will give examples of frontend codes
(which are in the main code folder). For editing the backend, such as adding
new geometries, sequences, or dynamics (such as flow), either check in my Honour's Thesis or email me. The design of the code
can handle any geometry, any sequence and any type of particle movement (not just diffusion) in any dimension, but not everything 
is implemented currently.

Finally, the UCL Camino group offers a monte carlo simulation code in their Camino library. However, this code
does not support multicompartment geometries with different diffusion coefficients when the permeability is nonzero
(see camino mailing list archives). They also use the high-field approximation always. These may not be of concern to a lot of people, but it is one extra feature that my code supports.

Compiling
---------

Terminology: the frontend code is the problem you want to solve. The backend code is the general
algorithms I have coded to solve general problems.

You will need an NVIDIA card that is FERMI(~2010) generation or newer.

1. Your frontend code file will typically look like "filename.cu" and will be located in the directory
you are compiling in (the main code directory).
2. type "cudaCompile filename" to compile the code in a cmd prompt. cudaCompile is a program I made so I wouldn't have to
remember all the junk in the cuda compiler command. If you do not remember how to use cudaCompile, you can always
type "cudaCompile" to see an example of its usage.
3. If you want to do it by entering the command, or you have extra include files which aren't in the current
directory, type this into a console:
"C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\bin\vcvars32.bat\" && "C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v5.5\bin\nvcc.exe\" --use-local-env --cl-version 2010 -ccbin "C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\bin\" --machine 32 -arch sm_20 -I "C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v5.5\include\" -I .\acquisition -I .\functors -I .\kernels -I .\params -I .\primitives -I .\substrate -I .\util -I .\blochdiff -o filename.exe filename.cu


Examples
--------

**********************************************************************

   Frontend Example 1 (GPU): Simulating Diffusion in a Single Object 

**********************************************************************

Geometry: Sphere
Sequence: PGSE

For an example with just the phase incrementing see magHF_spherNeuman.cu or spherePhase.cu
For an example with full magnetization evolution and concomitant gradients see magLF_spherNeuman.cu

Geometry: Parallel Planes
Sequence: Cos OGSE

For an example with just the phase incrementing see planePhase.cu

Geometry: Cylinder
Sequence: Cos OGSE

For an example with just the phase incrementing see cylPhase.cu

Geometry: Empty
Sequence: Cos OGSE

See Cos_Empty_Mel.cu

There are actually multiple (infinite) ways you can make an empty substrate in my framework.
The easiest would be probably to make a sphere/cylinder basis with zero radius (but be careful not to initialize
the particles in the sphere!). But you can also make your own empty basis class.

***********************************************************************************************************************

   Frontend Example 2 (GPU): Simulating Diffusion in a Multicompartment region with lattice periodicity and a cosOGSE
   (Gamma distribution of cylinders) 

***********************************************************************************************************************

See gammadistCos.cu. Gamma distribution parameters are from Matt Halls' paper http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=04797853
This code frontend has been heavily tested.

***********************************************************************************************************

   Frontend Example 3 (GPU): Complicated sequence with concomitant gradient terms 
   (UCN low-field gas experiment in a Big Cylinder) 

***********************************************************************************************************

See:
magLF_spherNeuman.cu
magLF_UCN.cu  (has a functor with an explicit 90 degree RF pulse!)
magLF_pgse.cu (has a functor with an explicit 90 and 180 degree RF pulse!)
cosSphere_LF.cu

for some examples 

