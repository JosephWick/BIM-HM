# Overview 
This repository features a hierarchical matrix implementation of a boundary integral method 2D SEAS simulation that includes viscoelastic effects. The ODE's for the time evolution of slip and strain are solved with Runge-Kutta 4/5 adaptive time stepping. 

The two main files are `BP1_visco_d.m` and `BP1_visco_hm.m`. Each solves the same physical problem but the former makes use of dense slip/strain kernels in a straighforward BIM implementation. The latter exchanges the dense kernels for hierachical matrices, which improve the matrix-vector product runtime for kernels with more than a few thousand elements. The physical problem is consistent with benchmark problem BP1-QD from the Southern California Earthquake Center Working Group for Advancing Simulations of Earthquakes and Aseismic Slip, which the addition of a viscoelastic region below the fault. The benchmark problems can be found [here](https://strike.scec.org/cvws/seas/benchmark_descriptions.html), and the modifications to this specific physical system are described below. 

The `include` directory includes modified versions of ode23 and ode45 that write data to disk to prevent the state vector from exceeding memory limits, as well as `hmmvp`, a tool for implementing hierachical matrices, and `ma2np`, a script that converts between matlab and python numpy arrays. More information on `hmmvp` can be found [here](https://github.com/ambrad/hmmvp), and `ma2np` can be found [here](https://github.com/joe-of-all-trades/mat2np)

# Physical Problem
The physical problem considered includes a fault -- with parameters matching the SCEC SEAS group's benchmark problem 1 -- and a viscoelastic zone. As can be seen in the below figure, the viscoelastic zone is 200 kilometers wide and stretches from just below the fault (a depth of 40km) to a depth of 200km.

![physProblem1 2](https://user-images.githubusercontent.com/39248450/177453087-9239cec4-7fed-4414-9b9e-bff9766425ef.png)
