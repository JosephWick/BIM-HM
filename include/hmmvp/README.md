# hmmvp-okada
Construct an H-matrix and compute matrix-vector products of the form B*x, B(rs,:)*x, and B(rs,cs)*x(cs).
***

## Installation
I used a distribution of CentOS linux, and unfortunately the Makefile included with `hmmvp` did not work for me. As a result, what follows is a less pretty way to get `hmmvp` installed and functional. I received a couple of warnings during compilation, but everything worked as intended. 

In the `include/hmmvp/` directory, first run `mkdir bin` and `mkdir lib`, then we can move onto compiling. We also need to compile `external/dc3omp.f` for linking with `hmmvp`. To compile the fortran code run the following in the `include/hmmvp/external/` directory:
```
gfortran -c dc3omp.f
```

Now we can go ahead and link everything into an hmmvp binary. In the `include/hmmvp` directory, run the following pieces of code. 

```
g++ -O3 -fopenmp -DUTIL_OMP -DFORTRAN_INT_4 -I . -c src/Hd.cpp -o src/Hd.o
g++ -O3 -fopenmp -DUTIL_OMP -DFORTRAN_INT_4 -I . -c src/Compress.cpp -o src/Compress.o
g++ -O3 -fopenmp -DUTIL_OMP -DFORTRAN_INT_4 -I . -c src/Hmat.cpp -o src/Hmat.o
g++ -O3 -fopenmp -DUTIL_OMP -DFORTRAN_INT_4 -I . -c src/HmatIo.cpp -o src/HmatIo.o
g++ -O3 -fopenmp -DUTIL_OMP -DFORTRAN_INT_4 -I . -c src/KeyValueFile.cpp -o src/KeyValueFile.o
g++ -O3 -fopenmp -DUTIL_OMP -DFORTRAN_INT_4 -I . -c src/CodeAnalysis.cpp -o src/CodeAnalysis.o
g++ -O3 -fopenmp -DUTIL_OMP -DFORTRAN_INT_4 -I . -c src/Mpi.cpp -o src/Mpi.o
g++ -O3 -fopenmp -DUTIL_OMP -DFORTRAN_INT_4 -I . -c src/CHmat.cpp -o src/CHmat.o
g++ -O3 -fopenmp -DUTIL_OMP -DFORTRAN_INT_4 -I . -c src/SFHmat.cpp -o src/SFHmat.o
ar rucs lib/libhmmvp_omp.a src/Hd.o src/Compress.o src/Hmat.o src/HmatIo.o src/KeyValueFile.o src/CodeAnalysis.o src/Mpi.o src/CHmat.o src/SFHmat.o
```
and 
```
gfortran src/hmmvpbuild.cpp external/dc3omp.o src/Hd.o src/Compress.o src/Hmat.o src/HmatIo.o src/KeyValueFile.o src/CodeAnalysis.o src/Mpi.o src/CHmat.o src/SFHmat.o -lstdc++ -fopenmp -llapack -lblas -o bin/hmmvpbuild_omp
```
If you want to get rid of the .o files, `make clean` should do the trick
