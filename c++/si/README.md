# Les Houches course: C++ code

# Shift-invert code

## Prerequesites

```
cd /path/to/home
mkdir src
cd src
```

1. Download the tar-ball of an older version of the linear solver Strumpack (newer versions have a default setup not favorable to shift-invert computations, this will be fixed at some point)
```
wget https://github.com/pghysels/STRUMPACK/archive/v2.2.0.tar.gz strumpack-2.2.0.tar.gz
```

2. Install, configure and compile PETSc

2a. Installation
```
git clone git@bitbucket.org:petsc/petsc.git
```

2b. Configuration : You will need cmake and a recent compiler. Here I use intel compiler, but gcc should be fine as well.  You should use the specific tarball of strumpack that you have downloaded. Pick a name for the PETSc ARCH, here I use real as PETSc is by default configured to work with real numbers.

Here is something that works for me on my local cluster
```
cd petsc
module load intel/19.4.243
module load intelmpi/19.4.243
module load cmake
./configure -with-debugging=0 --COPTFLAGS=O3 --CXXOPTFLAGS=O3 --FOPTFLAGS=O3 --with-blas-lapack-dir=$MKLROOT --with-errorchecking=0 --with-debugging=0 --COPTFLAGS=O3 --CXXOPTFLAGS=O3 --FOPTFLAGS=O3 --download-scalapack --download-parmetis --download-metis --download-strumpack=../strumpack-2.2.0.tar.gz --with-openmp PETSC_ARCH=real
```

2c. Compilation. Use the command line provided at the end of configuration. In my case it was:
```
make PETSC_DIR=path-to-home/src/petsc PETSC_ARCH=real
```

3. Install, configure and compile SLEPc
```
cd ..
git clone git@bitbucket.org:slepc/slepc.git
cd slepc

```

4. Adapt the ```CMakelists.txt` in the current directory to your needs. (we will add some help on this later).

5. ```cmake``` then `make spin_si`


## Running the code

'`spin_si` admits options which can either be given in command line (e.g. 'spin_si -L 12', ) or in the file `si.options`. Here is an exemple:
```
-L 16
-Sz 0
-disorder 5.0
-seed 3

-eps_nev 50
-st_type sinvert
-st_pc_factor_mat_solver_type strumpack
-st_pc_type lu
#-mat_strumpack_verbose
-mat_strumpack_colperm 0
-measure_entanglement 1
```
Command lines options take precedence over the ones in this file.
