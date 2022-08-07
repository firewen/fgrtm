# fgrtm
A code which uses the generalized reflection and transmission method (GRTM) to compute the synthetic seismogram

1, Installation
Before the installation, you should comfirm that IFORT or GFORTRAN have been installed. 
If the synthetics with SAC format is needed, you also install SAC first.
You should in the main directory and run "make". If there are any error information and an executable "MainProgram" appears in "bin",
that means the compilation succeed.

2, Running
There are one example in the "run".
In "single", you can change the parameters in GZ.BJT.dat according to your problem.
Then run "../../bin/MainProgram GZ.BJT.dat" to calculate the seismogram.
