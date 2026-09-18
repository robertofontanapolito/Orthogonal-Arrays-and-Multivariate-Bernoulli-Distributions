This repository contains the code required to reproduce the results presented in the paper "Orthogonal Arrays 
of Strength 2 for the Study of the Distributions of Pairwise-Independent Bernoulli(1/2) Random Variables", 
to appear in "Algebraic Statistics".

* the SAS script written in SAS/IML (`alg stat paper rev 01.sas`);

* the "4ti2" folder, containing:

  * input files for $d=3,4,5,6$, corresponding to the matrices $H_d$ (`.mat` files);
  * output files for $d=3,4,5$, containing the extreme rays (`.ray` files);
  * the commands used to run 4ti2 (`4ti2 command.txt`);

* the "oas" folder, containing:

  * the counting vectors of the OAs, stored as SAS dataset files. The filenames follow the convention
    `cnt_2_d_t2_N`, where `d` denotes the number of factors and `N` the size of the OA.

			 