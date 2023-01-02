# directedPercolation
This program measure the common observables of a directed percolation model and calculate critical exponents that belong to its universality class

To compile, run
cc directedPercolation.c my_nrutil.c -lm -O3 -o directedPercolation.x

To execute the program, run
./directedPercolation.x <Lattice size> <number of time interations> <percolation probability> <number of independent runs> <random negative integer>

After the program has been executed successfully, we have several ROOT macros to make graphs of the results. These are
  dataCollapse.c,
  fitting.c,
  fittingN_a.c,
  fittingP_sur.c,
  fittingR_s.c,
  makeMultiGraphs.c
