cmake_minimum_required(VERSION 2.8)

set(USE_MKL 1)

set(MKLROOT       $ENV{MKLROOT})
set(PETSC_DIR     $ENV{PETSC_DIR})
set(SLEPC_DIR     $ENV{SLEPC_DIR})
set(PETSC_ARCH     $ENV{PETSC_ARCH})
set(SLEPC_ARCH     $ENV{SLEPC_ARCH})

set(MPI_HOME      $ENV{I_MPI_ROOT})

set(LAPACK_INCLUDE $ENV{MKLROOT}/include)
set(LAPACK_LIB_DIR $ENV{MKLROOT}/lib/intel64/)
set(LAPACK_LIB     mkl_lapack95_lp64 mkl_intel_lp64 mkl_intel_thread mkl_core iomp5 pthread)

set(CMAKE_CXX_COMPILER "mpiicpc")
set(CMAKE_CC_COMPILER "mpiicc")
set(CMAKE_CXX_FLAGS "-O3 -qopenmp -fopenmp -std=c++0x -std=c++11")
set(CMAKE_C_FLAGS  "-O3 -qopenmp -fopenmp -std=c1x")

# set(MPILIB ${MPI_LIB})

include_directories(
${PETSC_DIR}/include/ ${SLEPC_DIR}/include/
${PETSC_DIR}/${PETSC_ARCH}/include/
${SLEPC_DIR}/${PETSC_ARCH}/include/
${MPI_HOME}/include/
${MPI_INC_DIR}
${LAPACK_INCLUDE}
)


link_directories(
        ${PETSC_DIR}/${PETSC_ARCH}/lib/
        ${SLEPC_DIR}/${PETSC_ARCH}/lib/
        ${MPI_HOME}/lib/
        ${MPI_LIB_DIR}
        ${LAPACK_LIB_DIR}
    )

### LAPACK, MKL
if(${USE_MKL})
add_definitions(-DUSE_MKL)
message("Using MKL")
endif(${USE_MKL})
message(${LAPACK_LIB_DIR} " using lapack: ${LAPACK_LIB}")


### EXEC FILEs
add_executable(spin_krylov spin_krylov.cc)
target_link_libraries(spin_krylov ${LAPACK_LIB} slepc petsc)
