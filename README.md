# leeds-code
Home of the Leeds Spherical Dynamo simulation code.

## Documentation
Documentation can be found in the online [Manual](https://leeds-spherical-dynamo.github.io/manual/)

Software dependencies
---------------------

The Leeds spherical dynamo code depends on recent versions of the following libraries:
- fftw (http://www.fftw.org/);
- blas (http://www.netlib.org/blas/);
- lapack (http://www.netlib.org/lapack/);
- netcdf (http://www.unidata.ucar.edu/software/netcdf/);

Compiling the parallel version of the code requires an MPI implementation, for example:
- openmpi (http://www.open-mpi.org),
- mpich (http://www.mpich.org);

These packages and the corresponding development headers should be installed in
standard locations or in locations that can be picked up by your compiler and linker.

Compilation and installation
----------------------------
On initially downloading the repository, run ./install.sh, to copy
parallel.h, Makefile, and parameters.F90 into their required paths, and edit them
to suit your environment/requirements.

If you are updating or changing branches, you should not need
to do this, so your locally edited Makefile, parallel.h and parameters.F90
can remain in place.

Once Makefile and parallel.h have been edited for your environment,
compile with `make`, and install with `make install`. To clear 
intermediate build modules/files, use `make clean`.  

To run, copy a statefile.cdf.in to your run folder with main.out. Optionally,
you may need codBC.cdf.in, compBC.cdf.in, or codBC.txt.in, codS.txt.in,
compBC.txt.in and compS.txt.in. 

Once all required files are in the run directroy, either run with
`mpirun -n X ./main.out >> OUT &` or submit a job script to the 
hpc scheduler from the run directory.

