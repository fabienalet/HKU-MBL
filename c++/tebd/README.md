# Les Houches lecture: "Numerical simulations for MBL"

# TEBD C++ code

## Huge disclaimer

This code has not been thoroughfully tested, if at all. Use at your own risks!
In particular, it is based on an older version of ITensor, so it is likely not to work out of the box. Nevertheless, it might be useful to look at it to understand the logics of both Itensor and TEBD.

## Prerequesites

This simple TEBD code for the XXZ random-field chain heavily uses ITensor which needs to be installed.
For this we can simply follow the instructions here:
http://www.itensor.org/docs.cgi?page=install&vers=cppv3


Let's try to use the newer Itensor v3.

```
cd /path/to/home
mkdir src
cd src
git clone https://github.com/ITensor/ITensor itensor
cd itensor
cp options.mk.sample options.mk
```

Now edit the ```options.mk``` file according to what is best. Tip: if you use MKL, adapting this file as this:
```
BLAS_LAPACK_LIBFLAGS=-L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_rt -lmkl_core -liomp5 -lpthread
BLAS_LAPACK_INCLUDEFLAGS=-I$(MKLROOT)/include
```
could simplify things. Then
```
make
```

Now go to the current directory (where this README.md file is). Adapt the first line of the Makefile, e.g.
```
LIBRARY_DIR=/path/to/home/src/itensor
```
and
```
make
```
will create an executable `tebd_XXZ`.
