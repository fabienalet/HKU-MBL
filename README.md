# HKU-MBL

# Repository associated to the lectures "Many-body Localization: introduction, computational challenges and open questions"

[Hong Kong Computational and Theoretical Physics Study Group 2023](https://quantummc.xyz/study-group/)

March 17th -29th, 2023

Fabien Alet (fabien.alet@cnrs.fr)

# Important Remark

Both notes and C++ codes date from 2019 and have not been re-tested again. There might be some very minor adaptations in case you use more recent version of Python or PETSc - SLEPc libraries. If you find something does not work anymore, please send an email, or better, do a PR with the correction. Thanks!


# Tutorials

There are 5 notebooks on the following topics (in order of complexity):

- Basic illustration of the Lanczos algorithm
- Static properties of the XXZ random field chain
- Dynamical properties of the XXZ random field chain
- Floquet Quantum circuits and time crystals
- Scars in the PXP and related models

These tutorials and the associated python codes have been written by Nicolas Macé, some time ago. Many thanks to him. Please drop him an email [  n [dot] mace [at] protonmail [dot] com ] if you like them and definitively contact him if you use the Python libraries in ```lib``` for your research!

See ```Setup.md``` for the setup for the Python Notebooks

# C++ codes

There are also three C++ codes that you should be able to use to perform large scale simulations of MBL problems (all for the random-field XXZ spin chain model, but pretty easy to adapt)

- Shift-invert code to get interior eigenpairs
- Krylov time evolution code
- Basic TEBD code

# Notes

There is a draft of lectures notes,  prepared for schools in Les Houches in 2019, and Cargese in 2022. They have been slightly refreshed and adapted for the study group. They will be further reshaped and updated in summer 2023. Stay tuned!
