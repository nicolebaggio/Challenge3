
# Challenge 3

In this code, I solve the Laplace equation modelling the heat diffusion over a square domain with a prescribed temperature on the whole boundary, by using the Jacobi iteration method.
In particular, I implement a parallel solver for the above problem for a generic f and solve the problem with f(x) = 8π2 sin(2πx) sin(2πy),whose exact solution is u(x, y) = sin(2πx) sin(2πy).





## Description of the code
The code only have the main.cpp file and the Makefile.
The main.cpp file has different utility function, needed for instance to evaluate L2-norm or to check the convergence.
The code asks to the user to set the dimension of the matrix (so the number of space-steps the method will performs) and the number of parallel tasks.
 
## OpenMP
The code also have proper OpenMP directive to further parallelize the local computations.


## Run the code and visualise results
The code have the file scalability.sh that runs a small scalability test (with 1,2, 4 cores ) and collect the results.
Moreover, there it is a file (plot_results.py) that plot the results.
To run tests, run the following commands

```bash
  mpicxx -fopenmp -o hybrid_parallel_program main.cpp
```

```bash
  mkdir data
```
```bash
  ./scalability_test.sh
```
```bash
  python plot_results.py

```

