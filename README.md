GMGPolar
=======

GMGPolar is a performant geometric multigrid solver using implicit extrapolation to raise the convergence order. It is based on meshes in tensor- or product-format. GMGPolar's focus applications are geometries that can be described by polar or curvilinear coordinates for which suited smoothing procedures have been developed.

If using GMGPolar, please cite:

M. J. Kühn, C. Kruse, U. Rüde. Implicitly extrapolated geometric multigrid on disk-like domains for the gyrokinetic Poisson equation from fusion plasma applications. Journal of Scientific Computing, 91 (28). Springer (2022). Link: https://link.springer.com/article/10.1007/s10915-022-01802-1

Tested plateforms
-----------------

Working
=======

* Linux x86_64 with GNU 9.3.0  compilers.    

Obtaining the source code
-------------------------

The GmgPolar Solver does not require any external libraries.
It is possible to link the code with the sparse direct solver ``MUMPS``.

* ``MUMPS`` is optionnal. To use it, compile the code with option -DUSE_MUMPS. It is 
  recommended to use the latest version (currently 5.4.1) but any version ulterior 
  to 5.1.0 should be okay. MUMPS is available freely on demand on the MUMPS consortium 
  website "mumps-solver.org".
	
The installation can be done by typing the following commands in your terminal

    # download the latest stable version
    # it will create a directory named GMGPolar

    git clone https://github.com/mknaranja/GMGPolar

Now that everything is ready, we can compile the solver.
Edit the file ``Makefile.in`` so that it reflects your configuration (path to libraries, file 
names, etc).


Building the library
--------------------
          
The build process is done using CMake:

    # Create build directory
    mkdir -p build && cd build
    # Configure
    cmake ..
    # Build
    cmake --build .

Currently, the default build process only supports gnu compiler although Intel compiler
has been successfully tested for some configurations.

Running GmgPolar
------------

You can run the solver without having to write a code (as we do in the next section). After building 
the library, a binary is created called ``./build/gmgpolar_simulation``, it takes parameters directly from command-line.

   
    # To try GmgPolar on a small problem size, without having to write any code,
    # ./build/gmgpolar_simulation uses default parameters with a grid 49 x 64.

    ./build/gmgpolar_simulation

    # For more details on the available parameters, see the main.cpp source code.
    # You can control the number of OpenMP threads used by changing the environment variable.
    # Note that only MUMPS is parallelized at the moment.

    export OMP_NUM_THREADS=4
  

Execution an example
-------------------------------------------------

Once the library is built, you can run the examples:

    # the verbose option defines the extent of the output

    ./build/gmgpolar_simulation --verbose 2

    # the option --debug 1 turns on internal debugging and compares the results of the C++ code 
    # with the results from the previous matlab implementation.
   
    ./build/gmgpolar_simulation --debug 1


Issue tracker
-------------
If you find any bug, didn't understand a step in the documentation, or if you
have a feature request, submit your issue on our
`Issue Tracker`: https://github.com/mknaranja/GMGPolar/issues
by giving:

- reproducible parameters
- computing environment (compiler, etc.)


Release Notes
-------------
* GmgPolar 1.0
1) Working multigrid cycle
2) In-house solver and possibility to link with MUMPS for the smoothing and coarse grid solution
3) Extrapolation strategies:
   
	a. No extrapolation (--extrapolation 0)

	b. Default implicit extrapolation (--extrapolation 1)

	c. Non-default implicit extrapolation with smoothing of all nodes on the finest level [experimental, use with care, convergence cannot be observed with residual] (--extrapolation 2)
6) Optimization of apply_A / build_rhs / apply_prolongation / build_Asc / apply_Asc_ortho
