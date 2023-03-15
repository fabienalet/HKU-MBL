
# Shift-invert code : timings

Using the following `si.options` file,

```
-Sz 0
-disorder 5.0
-seed 3

-st_type sinvert
-st_pc_factor_mat_solver_type strumpack
-st_pc_type lu
#-mat_strumpack_verbose
-mat_strumpack_colperm 0
-measure_entanglement 0
```

here are typical timings (run time in seconds) on a fairly recent CPU (with up to 36 MPI processes), using 1 OpenMP thread. This benchmark is from Summer 2019, you can probably do better by now!

| Number of MPI processes | 1 | 2 | 4 | 8 | 16 | 36 |
| --- | --- | --- | --- | --- | --- | --- |
| L=16 10 eigenv.   | 5s  |   |   |   |   |   |
| L=16 100 eigenv.   |  8s |   |   |   |   |   |
| L=18 10 eigenv.   |  25s | 17s  | 11s  | 8s  | 5s  | 4s  |
| L=18 100 eigenv.  |  50s | 33s  | 21s  | 13s  | 12s  | 7s  |
| L=20 10 eigenv.   | 297s  | 176s  | 141s  | 75s  | 50s  |  32s |
|  L=20 100 eigenv. |   | 299s  | 226s  | 120s  |  87s |  57s |


# Typical command lines :
```
./spin_si -L 16  -epsilon 0.5 -eps_nev 10
```
or
```
mpirun -np 36 ./spin_si -L 16  -epsilon 0.5 -eps_nev 10
```

would result in:
```
# L = 16 number of states=12870
# field= { -4.2927511954894335133 3.3994904246836625816 -3.7867142067036945186 0.69311325790087607857 -0.62938059705089077767 -4.8125198951543000092 -4.5936926241834061813 -2.5211169821972894134 -4.0644850440987685403 1.9482365726512149706 -3.5462985256922516797 -0.46827700348376932737 -2.8442298745523579839 -1.460949270337617012 -0.073640787330793600063 4.1330145289908699624 }
#Hamiltonian matrix assembly done.
Emin of H = -20.696176709417830608
Emax of H = 18.924395679713978069
#Processing target -0.88589051485192626956 [ epsilon=0.5 ] ... #Solved done.
E_0 = -0.88560479259787527884
E_1 = -0.88670801254207443076
E_2 = -0.88707821047332136022
E_3 = -0.88420151579507566542
E_4 = -0.88407118755815106148
E_5 = -0.88403356376629416147
E_6 = -0.88402778349104360611
E_7 = -0.88801326436184702118
E_8 = -0.88362070542304360643
E_9 = -0.88930174549111429272
E_10 = -0.88993689188449220939
E_11 = -0.89023919859335498561
E_12 = -0.88050030486092989435
*** Gap ratios
0.4759638282050380198
0.49294194455072015693
0.7257024315578480822
0.39591079807243645394
0.33556131140800132817
0.78617414753679970563
0.092874218874417380221
0.28868488322040081195
0.15363351127765934034
0.014199426854299014578
0.13045699098459692666
< r > = 0.35382759023111060026
```
