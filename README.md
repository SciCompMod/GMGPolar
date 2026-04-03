[![GMGPolarCI](https://github.com/SciCompMod/GMGPolar/actions/workflows/main.yml/badge.svg?branch=main)](https://github.com/SciCompMod/GMGPolar/actions/workflows/main.yml)
[![codecov](https://codecov.io/gh/SciCompMod/GMGPolar/graph/badge.svg?token=D0IVLUW51J)](https://codecov.io/gh/SciCompMod/GMGPolar)

![gmgpolar_logo](gmgpolar_small.png)

GMGPolar is a performant geometric multigrid solver using implicit extrapolation to raise the convergence order. It is based on meshes in tensor- or product-format. GMGPolar's focus applications are geometries that can be described by polar or curvilinear coordinates for which suited smoothing procedures have been developed.

If using GMGPolar, please cite:

M. J. Kühn, C. Kruse, U. Rüde. Implicitly extrapolated geometric multigrid on disk-like domains for the gyrokinetic Poisson equation from fusion plasma applications. Journal of Scientific Computing, 91 (28). Springer (2022). Link: https://link.springer.com/article/10.1007/s10915-022-01802-1

## Obtaining the source code

The GMGPolar Solver can run with or without the sparse direct solver `MUMPS`, though using MUMPS is recommended for optimal performance. This guide provides instructions on obtaining the code and installing the necessary dependencies.

## Clone the GMGPolar Repository

To begin, download the latest stable version of GMGPolar by running the following commands in your terminal:

```bash
# Clone the repository. This will create a directory named GMGPolar.
git clone https://github.com/SciCompMod/GMGPolar
```

## Configuring the Solver

After cloning the repository, you can configure the solver for your system. Configuration is mainly done by enabling options in `CMakeLists.txt` and setting paths to external libraries via environment variables for external dependencies (such as `MUMPS_DIR`, `METIS_DIR`, or `LIKWID_DIR`).

## Installing MUMPS using Spack

We highly recommend using Spack to manage and install external dependencies such as MUMPS. The following steps outline the process for installing MUMPS and related tools.

## Step 1: Install Spack

To install and set up Spack, execute the following commands in your terminal:

```bash
# Clone the Spack repository
git clone https://github.com/spack/spack.git

# Add Spack to your environment by sourcing its setup script
echo ". $HOME/spack/share/spack/setup-env.sh" >> ~/.bashrc

# Refresh your terminal or source your .bashrc
source ~/.bashrc
```

## Step 2: Install MUMPS

With Spack set up, you can now install MUMPS. The following command installs version 5.5.1 of MUMPS with specific options that are recommended for GMGPolar:

```bash
spack install mumps@5.5.1 ~blr_mt ~complex +double +float ~incfort ~int64 +metis ~mpi +openmp ~parmetis ~ptscotch ~scotch +shared
```

### Note on AVX / AVX-512 Compatibility

If your system does not support AVX or AVX-512 instructions (e.g., on AMD processors), install MUMPS with the following command:

```bash
spack install mumps@5.5.1 target=x86_64 ~blr_mt ~complex +double +float ~incfort ~int64 +metis ~mpi +openmp ~parmetis ~ptscotch ~scotch +shared
```

## Step 3: Configure CMake for GMGPolar

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

## Optional: Measuring performance

We use `Likwid` for performance monitoring. You can install it using Spack as well:

**Install Likwid (Performance Monitoring Tool)**:

```bash
spack install likwid
```

## Running GMGPolar

You can run the solver without having to write a code (as we do in the next section). After building
the library, a binary is created called `./build/gmgpolar`, it takes parameters directly from command-line.

    # To try GMGPolar on a small problem size, without having to write any code,
    # ./build/gmgpolar uses default parameters with a grid 33 x 64.

    ./build/gmgpolar

    # For more details on the available parameters, see the scripts/tutorial/run.sh script.

## Issue tracker

If you find any bug, didn't understand a step in the documentation, or if you
have a feature request, submit your issue on our
`Issue Tracker`: https://github.com/mknaranja/GMGPolar/issues
by giving:

- reproducible parameters
- computing environment (compiler, etc.)

## Release Notes

### GMGPolar 1.0.0

1. **Multigrid**
    - Implemented a fully functional multigrid cycle
3. **Direct Solver**
    - Added an in-house direct solver
    - Enabled optional integration with MUMPS for smoothing and coarse-grid solves
5. **Extrapolation strategies**
   - No extrapolation
   - Default implicit extrapolation
   - Non-default implicit extrapolation with smoothing of all nodes on the finest level [experimental, use with care, convergence cannot be observed with residual]
6. **Optimization**
   - Improved performance of key routines: apply_A / build_rhs / apply_prolongation / build_Asc / apply_Asc_ortho

### GMGPolar 2.0.0

1. **Enhancements and New Class Layout**

- **Linear Algebra:**
  - Introduced custom Vector and SparseMatrix classes.
  - Added a (cyclic) Tridiagonal Solver for improved performance and usability.
- **Input Functions:**
  - Separated into distinct components: DomainGeometry, BoundaryConditions, SourceTerm, etc.
- **Polar Grid:**
  - Indexing is now based on circle/radial smoother.
- **Residual:**
  - Improved the residual calculation by addressing the unintuitive behavior that previously applied only to the interior part of the matrix.
- **Direct Solver:**
  - Fixed a bug where boundary values were not treated correctly.
  - Built matrices to be symmetric, reducing factorization time.
- **Smoother:**
  - Separated into extrapolated and standard smoothers.
  - Replaced the LU-Decomposition algorithm with the Thomas algorithm for improved efficiency.

2. **New Features**
    - Introduced W- and F cycles for enhanced solving capabilities.
    - Added FMG (Full Multigrid) to obtain improved starting solutions.
    - Implemented advanced caching behavior options for the "A-Give" implementation strategy.
    - Added a new strategy named 'A-Take,' which is appropriate for cases where memory is less of a constraint, resulting in faster execution times.
    - Comprehensive Unit Tests: Integrated Google Unit Tests for all classes, ensuring robust and reliable functionality across the codebase.

3. **Performance Improvements**
    - Removed the task-based approach, which did not scale well with increasing parallelization.
    - Reduced maximum usage by 61.5% by constructing symmetric matrices and utilizing the tridiagonal structure of smoother matrices.

4. **Updated Features**
    - Added a new LU decomposition solver, allowing users to choose between MUMPS and the in-house solver for greater flexibility and performance.

### GMGPolar 2.0.1

1. **Minor Changes**
    - Correction of broken code coverage
    - Improved output
    - Changed position for anisotropic refinement with ZoniShifted

### GMGPolar 2.1.0

1. **Solver & Performance**
    - Replaced the custom LU decomposition solver with a faster implementation.
    - Fixed an error in the FMG method, reducing iterations when using a small number of multigrid levels.

2. **New Functionality**
    - Added angular dependence in profile coefficients.
    - Enabled multiple solves with the same GMGPolar instance for varying source terms and boundary conditions.

3. **Restructuring & Code Quality**
    - Refactored GMGPolar by moving command-line parsing into a dedicated ConfigParser.
    - Redesigned the solve and setup functions, and reorganized the header file with a clearer structure and enhanced documentation.

4. **Testing & Usability**
    - Added tests for LU solver, convergence order, and solver validation across multiple settings.
    - Improved verbose output formatting for clearer settings and runtime information.

### GMGPolar 2.2.0

1. **Preconditioned Conjugate Gradient (PCG)**
   - Added PCG solver, allowing GMGPolar to be used as a preconditioner for CG instead of a standalone iterative solver.
   - When solving the extrapolated problem, PCG converges in up to 4x fewer iterations and runs up to 2x faster end-to-end.
   - Addtitional memory overhead is minimal by aliasing PCG work vectors onto existing storage.

2. **Kokkos Integration**
   - Switched to Kokkos-backed vectors and a batched tridiagonal solver with Kokkos parallelization.
   - Templated multigrid operators to support future GPU execution spaces, removing the previous polymorphic design.

3. **Design and API Cleanup**
   - Input functions now use C++ concepts instead of polymorphism.
   - Simplified interpolation class to interpolate directly between grids instead of levels.
   - Encapsulated MUMPS solver in its own class; improved DirectSolver naming consistency.
   - Replaced macro-heavy patterns with standard functions for type safety.
   - Removed unused Point, MultiIndex classes, redundant LevelCache constructor, and thread reduction factor variable.

4. **New Features**
   - Added support for solves without a multigrid hierarchy.

5. **Bug Fixes**
   - Fixed MUMPS factorization failure when OpenMP multithreading is enabled in versions later than 5.5.1. Version 5.5.1 remains recommended.

6. **Testing**
   - Added formatting validation tests and automatic CI testing with MUMPS.
   - Added Google Tests for PCG convergence validation.
