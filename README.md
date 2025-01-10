GMGPolar
=======

GMGPolar is a performant geometric multigrid solver using implicit extrapolation to raise the convergence order. It is based on meshes in tensor- or product-format. GMGPolar's focus applications are geometries that can be described by polar or curvilinear coordinates for which suited smoothing procedures have been developed.

If using GMGPolar, please cite:

M. J. Kühn, C. Kruse, U. Rüde. Implicitly extrapolated geometric multigrid on disk-like domains for the gyrokinetic Poisson equation from fusion plasma applications. Journal of Scientific Computing, 91 (28). Springer (2022). Link: https://link.springer.com/article/10.1007/s10915-022-01802-1

## Obtaining the source code

The GMGPolar Solver requires external libraries, specifically the sparse direct solver ``MUMPS``. This guide provides instructions to obtain the code and install the necessary dependencies.

## Clone the GMGPolar Repository

To begin, download the latest stable version of GMGPolar by running the following commands in your terminal:

    # Clone the repository. This will create a directory named GMGPolar.
    git clone https://github.com/mknaranja/GMGPolar

## Configuring the Solver

After cloning the repository, you'll need to configure the solver for your system. Edit the ``CMakeLists.txt`` file to reflect your system's configuration (e.g., paths to libraries, file names, etc.).

## Installing MUMPS using Spack

We highly recommend using Spack to manage and install external dependencies such as MUMPS. The following steps outline the process for installing MUMPS and related tools.

## Step 1: Install Spack

To install and set up Spack, execute the following commands in your terminal:

    # Clone the Spack repository
    git clone https://github.com/spack/spack.git

    # Add Spack to your environment by sourcing its setup script
    echo ". $HOME/spack/share/spack/setup-env.sh" >> ~/.bashrc

    # Refresh your terminal or source your .bashrc
    source ~/.bashrc

## Step 2: Install MUMPS

With Spack set up, you can now install MUMPS. The following command installs version 5.5.1 of MUMPS with specific options that are recommended for GMGPolar:

    spack install mumps@5.5.1
        blr_mt=false
        complex=false
        double=true
        float=true
        incfort=false
        int64=false
        metis=true
        mpi=false
        openmp=true
        parmetis=false
        ptscotch=false
        scotch=false
        shared=true

## Step 3: Install Required Dependencies

MUMPS relies on other packages like `Metis` for matrix ordering and `Likwid` for performance monitoring. You can install these dependencies using Spack as well:

1. **Install Metis (Fill-Reducing Matrix Ordering Library)**:
    ```bash
    spack install metis
    ```

2. **Install Likwid (Performance Monitoring Tool)**:
    ```bash
    spack install likwid
    ```

## Step 4: Configure CMake for GMGPolar

After installing MUMPS and other dependencies, ensure that the paths to the libraries are correctly set in the CMakeLists.txt file.

## Final Step: Compiling the GMGPolar Solver

Once everything is configured, compile the solver by running the following commands:

```bash
cd scripts
./compile.sh [Debug|Release]
```

After executing ./compile.sh [Debug|Release], the script will compile the solver using the specified build type. You can also run ./compile.sh without any arguments afterward, and it will automatically use the last configured build type.

Currently, the default build process only supports gnu compiler although Intel compiler
has been successfully tested for some configurations.

## Running GMGPolar

You can run the solver without having to write a code (as we do in the next section). After building 
the library, a binary is created called ``./build/gmgpolar``, it takes parameters directly from command-line.

    # To try GMGPolar on a small problem size, without having to write any code,
    # ./build/gmgpolar uses default parameters with a grid 33 x 64.

    ./build/gmgpolar

    # For more details on the available parameters, see the scripts/run.sh script.
  
## Issue tracker

If you find any bug, didn't understand a step in the documentation, or if you
have a feature request, submit your issue on our
`Issue Tracker`: https://github.com/mknaranja/GMGPolar/issues
by giving:

- reproducible parameters
- computing environment (compiler, etc.)

## Release Notes

## GMGPolar 1.0.0
1) Working multigrid cycle
2) In-house solver and possibility to link with MUMPS for the smoothing and coarse grid solution
3) Extrapolation strategies:
   
	a. No extrapolation (--extrapolation 0)

	b. Default implicit extrapolation (--extrapolation 1)

	c. Non-default implicit extrapolation with smoothing of all nodes on the finest level [experimental, use with care, convergence cannot be observed with residual] (--extrapolation 2)
6) Optimization of apply_A / build_rhs / apply_prolongation / build_Asc / apply_Asc_ortho


## GMGPolar 2.0.0

### Enhancements

- **Refactored Class Layout:** 
  - **Linear Algebra Enhancements:**
    - Introduced custom **Vector** and **SparseMatrix** classes.
    - Added a (cyclic) **Tridiagonal Solver** for improved performance and usability.
  - **Input Function Enhancements:**
    - Separated into distinct components: **DomainGeometry**, **BoundaryConditions**, **SourceTerm**, etc.
  - **Polar Grid:**
    - Indexing is now based on circle/radial smoother.
  - **Residual**
    - Improved the residual calculation by addressing the unintuitive behavior that previously applied only to the interior part of the matrix.
  - **Direct Solver:**
    - Fixed a bug where boundary values were not treated correctly.
    - Built matrices to be symmetric, reducing factorization time.
  - **Smoother:**
    - Separated into extrapolated and standard smoothers.
    - Replaced the LU-Decomposition algorithm with the Thomas algorithm for improved efficiency.
  
### New Features

- Introduced **W** and **F** cycles for enhanced solving capabilities.
- Added **FMG** (Full Multigrid) to obtain improved starting solutions.
- Implemented advanced caching behavior options for the "Give" implementation strategy.
- Added a faster strategy named 'Take,' which is appropriate for cases where memory is less of a constraint, resulting in an 80% increase in memory usage.
- **Comprehensive Unit Tests:** Integrated Google Unit Tests for all classes, ensuring robust and reliable functionality across the codebase.

### Performance Improvements

- Removed the task-based approach, which did not scale well with increasing parallelization.
- Reduced maximum usage by 61.5% by constructing symmetric matrices and utilizing the tridiagonal structure of smoother matrices.

### Removed Features

- Discontinued the in-house solver, now replaced by MUMPS for improved functionality.