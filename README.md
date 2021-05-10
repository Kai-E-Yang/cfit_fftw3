# cfit
## What to do
Calculating the non-linear force foree field with current iteration based on the Grad-Rubin algorithm. This is used in the solar physics. Using the bottom magnetic field information from the solar observations to reconstruct the 3D magnetic field in the solar corona. If you use it, please cite the following papers,
```
@ARTICLE{2006SoPh..238...29W,
       author = {{Wheatland}, M.~S.},
        title = "{A Fast Current-Field Iteration Method for Calculating Nonlinear Force-Free Fields}",
      journal = {\solphys},
         year = 2006,
        month = oct,
       volume = {238},
       number = {1},
        pages = {29-39},
          doi = {10.1007/s11207-006-0232-0},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2006SoPh..238...29W},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

@ARTICLE{2007SoPh..245..251W,
       author = {{Wheatland}, M.~S.},
        title = "{Calculating and Testing Nonlinear Force-Free Fields}",
      journal = {\solphys},
     keywords = {Active regions: magnetic fields, Active regions: models, Corona: models, Magnetic fields: corona, Magnetic fields: models},
         year = 2007,
        month = oct,
       volume = {245},
       number = {2},
        pages = {251-262},
          doi = {10.1007/s11207-007-9054-y},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2007SoPh..245..251W},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```
## What in this program
- The file main.f03 is the main program.
- The file mod_io.f03 contains all the routines used for reading and writing.
- The file mod_param.f03 contains all the definitions of global parameters.
- The file mod_operator.f03 contains all the small parts of some calculations, like rk4, cross product ...
- The file mod_solver.f03 contains the main pathway for the calculation.
- The file mod_fftw.f03 contains the subroutines for FFTW, and the Poisson, Laplace solver based on FFTW.
- This program is used for calculating the squashing factor Q.
- The method is proposed by Wheatland 2006, 2007, SolPhy.

## History:
- The code was modified by Xue Yang (Jan 2008) to use the Intel MKL Fourier routines, which do not require padding to a power of two. 
- Versions 8-10: Padding removed. (If required, padding should be included in the BCs supplied to the code.) Progressive censoring removed. Mixed BCs option (ABC=0) removed. Periodic field line tracing implemented (flag PER). Self-consistency cycling implemented (flag NCYC). Closed top BC implemented (flag TBC). 
- Version 11: Uncertainties in alpha included, for use in self-consistency procedure.
- Version 12:
- Version 13: Different handling of uncertainties in alpha, and CONV flag introduced.
- Version 14: Removed AF, CONV flags, included use of periodicity for interpolation of points leaving the box.
- Version 15: Introduced least-squares fitting of alpha values in the volume (LSQ flag)
- Version 20: Code parallelised with mpi (Stuart Gilchrist 2012). 
-- Method of gathering results from MPI nodes has been improved. 
-- Dynamic scheduling is now used with the OpenMP tracing loop parallelisation. 
- Version 20.3 (changes):
-- Change to weightings in trapezoidal rule in the energy() routine in check.f. The new weightings correspond to Gaussian quadrature for the periodic volume.
- Version 21: Implemented flux-balanced top boundary. If TBC true and total bottom flux is nonzero, Bz at the top will be uniform, to balance out the bottom flux        
- Version 22.1: Can now accommodate iterations up to 9999
- Version 23: Can now use LSQ
- Version 26: 2021 May 07, changed by Kai E. Yang.
-- change style of the code.
-- used FFTW3 for the Fourier Transfer instead of Intel MKL FFT.
-- change the keywords' name.
-- change the input format of the keywords as the fortran namelist.
-- change the log field into CSV format.

## How to use it
1, set the directory of the FFTW3 lib, e.g.,
```
FFTW3_DIR=/opt/fftw3
```
The current version FFTW3.3.3 has been tested.

2, compile the program by using the makefile
```
>make -f Makefile
```

3, create 2d bottom boundary conditions for Bn, alpha, (optionally, sigma of the alpha), example codes in IDL and python have been provided in ./test1 and ./test2. Test1 provide the Sturrock linear force-free field as the test field. Test2 provide the Low & Lou non-linear force-free field as the test field.
```
> IDL
IDL> .r bcs.pro
```
or
```
>python bcs.py
```

4, make a soft link of the program at the target directory, like ./test1
```
cd ./test1
ln -s ../code/cfit
```

5, run the code with parameter file "par"
```
./cfit input.par
```

5, postprocessing codes has been proovided ./io, you can used them to read and visulize the results, examples have been provided in ./test1 and ./test2
```
>python read_cfit.py
```

## Parameter file
The format of the input.par file is the namelist in fortran, e.g.,

```
&FILENAME_PAR
    OUTFILENAME="./result/",
    BZ0NAME="bz0.dat",
    INDATAFORMAT="binary",
    ALPHANAME="alpha0.dat",
    ALPHAERRNAME="alpha0err.dat",
    RESTARTNAME="restartB.dat",
 /
&CAL_PAR
    NUMTHREADS=1,
    DIMX=128,
    DIMY=128,
    DIMZ=128,
    DELTA_S= 0.5d0,
    METHOD="fft",
    STEPMETHOD="rk2",
    NCYCLE=1,
    CHECK=F,
    RESTART=F,
    SAVNUM=1,
    FACTOR= 0.5d0,
    TOP_CLOSED=T,
    SPK_FLAG=F,
    SPK=  1.0E-005,
    PERIODIC=T,
    AOUT=T,
    POLARITY=1,
    NLOOP=1,
    STARTCYC=1,
    STARTLOOP=1,
    ALPHA_ERROR=F,
    MIN_STEP=5.0E-002,
 /
```
The above are default settings for the keywords, which is used to control the calculation.
The keywords are classified as two groups, the group FILENAME_PAR contains all the keynames for the io file and data format. The group CAL_PAR contains all the keywords to control the calculations.
- FILENAME_PAR
-- OUTFILENAME, the directory of the output result.
-- BZ0NAME, the file name of BCs of the Bz value.
-- INDATAFORMAT, data format.
-- ALPHANAME, file name of BCs of the alpha value
-- ALPHAERRNAME, file name of the BCs of the sigma of alpha.
-- RESTARTNAME, name of the restart file that contains a 3D B field.
- CAL_PAR
-- NUMTHREADS, the number of threads used for the parallization.
-- DIMX, numbers in X dimension.
-- DIMY, numbers in Y dimension.
-- DIMZ, numbers in Z dimension.
-- DELTA_S, step length for the field line integral, in the unit of pixel.
-- METHOD, a string indicate the mathod for Poisson and Laplace eqaution.
-- STEPMETHOD, method for the field line integral.
-- NCYCLE, how many cycles for the self-consistent solution, if 1, indicate NO self-consistent solution.
-- CHECK, a logical flag to determine whether output the electric current, left_box (indicate the field line is closed or not), 3d alpha value assigned by the field line tracing, and the electric current associated magnetic field.
-- RESTART, if this parameter is set True, the program will restart from an existed magnetic field, stored in the file RestartName.
-- SAVNUM, the step number for saving a snapshot of the calculation. Actually, the default value is set as the same as the NCYCLE.
-- FACTOR, the factor to update the magnetic field, which makes the calculation stable, i.e., B_new=B_new*factor+B_old*(1-factor). It means the larger factor, the faster the magnetic field updated.
-- TOP_CLOSED, if it is true, the top boundary is chosen as a closed boundary, i.e., Bz at top = 0. Othervise the Bz will vanish at infinity.
-- SPK_FLAG, if true, the program will reduce the spikes in the Bx and By on the bottom.
-- SPK, the value of the reducing spikes.
-- PERIODIC, if true, the field lines will be traced periodically on the side boundaries.
-- AOUT, if true, the associated vector potential of the magnetic field will be written in OUTFILENAME.
-- POLARITY, the value of this parameter is 1 or -1, 1 means the initial cycle is the possitive solution, -1 means the initial cycls is negative solution.
-- NLOOP, the iteration number in each polarity cycle.
-- STARTCYC, if the program is restart from some existed field, the STARTCYC can be changed, which reset the numbers in the written results. It should be smaller than NCYCLE.
-- STARTLOOP, if the program is restart from some existed field, the STARTLOOP can be changed, which reset the numbers in the written results. It should be smaller than NLOOP.
-- ALPHA_ERROR, this logical flag indicate whether use the sigma values stored in ALPHAERRNAME
-- MIN_STEP, minimum step size of the field line tracing in the unit of pixel.

## Other Note:
This parallel version is based on FORTRAN OPENMP.

If one want to use N threads for the calculation, just change the value of parameter 'nthreads' in parameter file. If this parameter is defined as 0, then the max number of threads in the computer will be used as default. The default value in the program is 1.

The field line integral method is Runge-Kutta 4(5).

## TODO
- Impliment the MPI for the parallization.
