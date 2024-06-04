# Challenge 3

In this code, I solve the Laplace equation modeling heat diffusion over a square domain with a prescribed temperature on the whole boundary, using the Jacobi iteration method. Specifically, I implement a parallel solver for the above problem for a generic function f and solve the problem with f(x) = 8π2 sin(2πx) sin(2πy),whose exact solution is u(x, y) = sin(2πx) sin(2πy).

## Description of the Code

The code consists of the following main files: `main.cpp` ,`chrono.hpp`, `Makefile`, `scalability_test.sh`, `plot_results.py`.

- **main.cpp**: Contains the main implementation along with utility functions for tasks such as evaluating the L2-norm and checking convergence.
- **chrono.hpp**: Is an utility to evaluate the time used by the code to perform.
- **Makefile**: Provides a way to compile the code easily.
- **scalability_test.sh**: Runs a small scalability test (with 1,2, 4 cores ) and collect the results.
- **plot_results.py**: Plots the execution time and the L2-norm-error over the number of cores.


The user can set the dimension of the matrix (i.e., the number of space steps the method will perform) and the number of parallel tasks by modifying the file `scalability_test.sh`, by changing this part of the code:
```scalability_test.sh
   # Number of processors to test
   PROCS=(1 2 4)

   # Grid size
   N=64
```


## OpenMP

The code includes proper OpenMP directives to further parallelize local computations.

## Run the Code and Visualize Results

To run the tests, follow these commands:

1. **Compile the Program**:
```bash
   mpicxx -fopenmp -o hybrid_parallel_program main.cpp
```
2. **Create the Data directory**:
```bash
  mkdir data
```
3. **Run the Scalability Test**:
```bash
  ./scalability_test.sh
```
4. **Generate and view the Performance Plot**:
```bash
  python plot_results.py

```

## Alternative main.cpp
As an alternative, in the main.cpp you can find a piece of commented code that permits the user to set the size of the matrix and the number of processors interactively. In this way you can run the code with the follow commands:

```bash
   make
```
```bash
   ./main
```