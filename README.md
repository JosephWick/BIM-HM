# Overview 
This repository features a hierarchical matrix implementation of a boundary integral method 2D SEAS simulation that includes viscoelastic effects. The ODE's for the time evolution of slip and strain are solved with Runge-Kutta 4/5 adaptive time stepping. 

The two main files are `BP1_visco_d.m` and `BP1_visco_hm.m`. Each solves the same physical problem but the former makes use of dense slip/strain kernels in a straighforward BIM implementation. The latter exchanges the dense kernels for hierachical matrices, which improve the matrix-vector product runtime for kernels with more than a few thousand elements. The physical problem is consistent with benchmark problem BP1-QD from the Southern California Earthquake Center Working Group for Advancing Simulations of Earthquakes and Aseismic Slip, which the addition of a viscoelastic region below the fault. The benchmark problems can be found [here](https://strike.scec.org/cvws/seas/benchmark_descriptions.html). This specific problem differs from the standard benchmark one problem in that it includes a viscoelastic region that extend 200km wide, and begins right below the fault at a depth of 40km and extends down to a depth of 200km. 

The `include` directory houses modified versions of ode23 and ode45 that write data to disk to prevent the state vector from exceeding memory limits, as well as `hmmvp`, a tool for implementing hierachical matrices, and `ma2np`, a script that converts between matlab and python numpy arrays. More information on `hmmvp` can be found [here](https://github.com/ambrad/hmmvp), and `ma2np` can be found [here](https://github.com/joe-of-all-trades/mat2np)

# Code Specifics
## Usage
The dense script can be run in a single step, by calling it in a matlab console. The hierarchical version of the simulation is a three step process, requiring first a setup to be ran, then creation and compression of the hierarcical matrices, followed by running the actual simulation. Once kernels have been made, the simulation can be re-run without the first two steps so long as the .hm files remain available and no changes to the physical problem have been made. 

To run the setup, call `b = BP1_visco_hm("build")` in a matlab console. This create the starter files for hmmvp and store in memory the mesh information as the matlab variable `b`. 

Additionally, several terminal commands will be printed. These need to be ran in the same directory in a shell and will make use of hmmvp tools to create the compressed hierarchical matrices for later use. 

Once the build function is called and the hierarchical matrices are created, call `BP1_visco_hm("run", b)` in a matlab console. This will run the simulation according to the physical problem described by the build function. Run will also return references to the hierarchical matrices as a matlab structure. 

## HMMVP
### Installation
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

### Creating Hierarchical Matrices
Data is passed between matlab and hmmvp using key value files, denoted with the extension `.kvf`. Strings, numbers and arrays can all be passed through key value files, but for this project we only need to pass along arrays of floats. An example of writing to a key value file is below
```
% parameters for hmmvp
c.greens_fn = 'okadaS12';
c.write_hmat_filename = './tmp/BP1v_ff-s12';
c.kvf = [c.write_hmat_filename '.kvf'];

% parameters for the physical setup
c.Y = [faultX; faultY; faultZ'];
c.X = [faultX_c; faultY_c; faultZ_c];
c.L = Lf';
c.W = Wf';

% write .kvf file
kvf('Write', c.kvf, c, 4);
```

The parameters for hmmvp specify the file that defines the function used in creating the kernel, as well as filenames. Physical parameters are the information of fault mesh blocks that are passed along to hmmvp when it builds the kernel. 

Creating a custom function for hmmvp to use during kernel creation is rather straightforward. One must add their new function to the `hmmvpbuild.cpp` file, as well as create a class file for it. There are several custom functions in the `include/hmmvp/src/` directory. 

One item of note is the size of the hierarchical matrices. In the class declaration of each custom function you will see the line:
`virtual Hd* ComputeHd (double eta) { return NewHd(_z, _x, NULL, eta); }`
`_z` and `_x` are parameters passed into the function by way of key value files. They are also what is used to define the shape of the kernel. Hmmvp takes the length of the first array contained in `_x` for the number of columns in the kernel, and the length of the first array contained in `_z` for the number of rows. 

In the case where you are considering self interaction (e.g. the fault interacting with itself) the kernels are square, and creating properly sized kernels is generally trivial. In the case of rectangular kernels (such as the fault feeling effects of the shear zone), I found it easiest to create a matrix to pass into hmmvp that has no use but for sizing. In the code this is depicted as `_z` or just `Z` in matlab. 
