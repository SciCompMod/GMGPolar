GMGPolar
=======

GMGPolar is a performant geometric multigrid solver using implicit extrapolation to raise the convergence order. It is based on meshes in tensor- or product-format. GMGPolar's focus applications are geometries that can be described by polar or curvilinear coordinates for which suited smoothing procedures have been developed.

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
          
The build process is done using ``make``:

    # run make
    make

The code is automatically built with `build_intel` or `build_gnu`, depending on your compiler.

Running GmgPolar
------------

You can run the solver without having to write a code (as we do in the next section). After building 
the library, a binary is created called ``./build/main``, it takes parameters directly from command-line.

   
    # To try GmgPolar on a small problem size, without having to write any code,
    # ./build/main uses default parameters with a grid 49 x 64.

    ./build_gnu/main

    # For more details on the available parameters, see the main.cpp source code.
    # You can control the number of OpenMP threads used by changing the environment variable.
    # Note that only MUMPS is parallelized at the moment.

    export OMP_NUM_THREADS=4
  

Building an example (to call ABCD from C++ or C)
-------------------------------------------------

Once the library is built, you can run the debug examples (either C++ or C):

    # the option --debug 1 turns on debugging and compares the results of the C++ code 
    # with the results from the previous matlab implementation.
   
    ./build/main --debug 1


Issue tracker
-------------
If you find any bug, didn't understand a step in the documentation, or if you
have a feature request, submit your issue on our
`Issue Tracker <https://github.com/mknaranja/GMGPolar/issues>`
by giving:

- reproducible parameters
- computing environment (compiler, etc.)


Release Notes
-------------
* GmgPolar 0.9
1) Working multigrid cycle
2) In-house solver and possibility to link with MUMPS for the smoothing and coarse grid solution
3) Extrapolation strategies:
    a. No extrapolation (--extrapolation 0)
    b. Implicit extrapolation with smoothing of fine nodes only (--extrapolation 1)
    c. Extrapolation with smoothing of all nodes (--extrapolation 2)
4) Optimization of apply_A / build_rhs / apply_prolongation / build_Asc / apply_Asc_ortho
