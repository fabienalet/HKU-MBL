# Les Houches lecture: "Numerical simulations for MBL"

# Krylov C++ code

**Disclaimer** : this code hasn't be thoroughfully tested, do not take results for granted.


1a. Configure PETSc with complex numbers.
If you have already installed Petsc, you just need to create a different ARCH. (we can simplify the configuration line as there is no need of external linear solvers here for the krylov code).

If you haven't, please look at README.md from the shift-invert code.

```
cd src/petsc
module load intel/19.4.243
module load intelmpi/19.4.243
./configure -with-debugging=0 --COPTFLAGS=O3 --CXXOPTFLAGS=O3 --FOPTFLAGS=O3 --with-blas-lapack-dir=$MKLROOT --with-errorchecking=0 --with-debugging=0 --COPTFLAGS=O3 --CXXOPTFLAGS=O3 --FOPTFLAGS=O3 --with-scalar-type=complex PETSC_ARCH=complex
make PETSC_DIR=path-to-home/src/petsc PETSC_ARCH=complex
```


1b. Compilation. Use the command line provided at the end of configuration. In my case it was:
```
make PETSC_DIR=path-to-home/src/petsc PETSC_ARCH=complex
```

2. Install, configure and compile SLEPc
```
cd ..
git clone git@bitbucket.org:slepc/slepc.git
cd slepc

```

3. Adapt the ```CMakelists.txt``` in the current directory to your needs. (we will add some help on this later).

4. ```cmake``` then `make spin_krylov`


## Running the code

`spin_krylov` admits options which can either be given in command line (e.g. 'spin_si -L 12', ) or in the file `krylov.options`. Here is an example:
```
-L 16
-Sz 0
-disorder 5.0
-seed 3
-num_times 100
# -loggrid 1
-Tmin 0
-Tmax 100
-cdw_start 1
#-product_state_start 1
#-num_product_states 20

-measure_imbalance 1
-measure_entanglement 1
#-measure_local 1
```

A detailed description of options will come later.  Command lines options take precedence over the ones in this file.
