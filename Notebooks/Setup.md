# Setup for the Python Notebooks

1. Install anaconda or miniconda
2. Create a new environment
```
conda create -n cargese intelpython3_core python=3
```
and activate it.

3. Install jupyter
```
conda install -c conda-forge jupyterlab
```
4. Install MKL
```
conda install -c intel intelpython3_full
```
5. Install quspin with OMP support
```
conda install -c weinbe58 omp quspin
```
