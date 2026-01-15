#pragma once

/* ---------------------------------- */
/* GMGPolar - Enumeration Definitions */
/* ---------------------------------- */

enum class StencilDistributionMethod
{
    CPU_TAKE = 0,
    CPU_GIVE = 1
};

/* Multigrid Cycle Types */
enum class MultigridCycleType
{
    V_CYCLE = 0,
    W_CYCLE = 1,
    F_CYCLE = 2
};

/* Residual Norm Type */
enum class ResidualNormType
{
    EUCLIDEAN          = 0, // Corresponds to the L2 norm
    WEIGHTED_EUCLIDEAN = 1, // Scaled L2 norm (sqrt(u^T * u / DoFs))
    INFINITY_NORM      = 2 // Corresponds to the L∞ norm
};

enum class ExtrapolationType
{
    NONE                         = 0,
    IMPLICIT_EXTRAPOLATION       = 1,
    IMPLICIT_FULL_GRID_SMOOTHING = 2,
    COMBINED                     = 3,
};

/* Smoother Colors */
enum class SmootherColor
{
    Black = 0,
    White = 1,
};

/* -----------*/
/* Test Cases */
/* -----------*/

/* Geometry Types - domain_geometry */
enum class GeometryType
{
    CIRCULAR  = 0,
    SHAFRANOV = 1,
    CZARNY    = 2,
    CULHAM    = 3
};

/* Test Problem Types - exact_solution */
enum class ProblemType
{
    CARTESIAN_R2   = 0,
    CARTESIAN_R6   = 1,
    POLAR_R6       = 2,
    REFINED_RADIUS = 3
};

/* Alpha Coefficient Types - profile_coefficient alpha */
enum class AlphaCoeff
{
    POISSON       = 0,
    SONNENDRUCKER = 1,
    ZONI          = 2,
    ZONI_SHIFTED  = 3
};

/* Beta Coefficient Types - profile_coefficient beta */
enum class BetaCoeff
{
    ZERO          = 0,
    ALPHA_INVERSE = 1
};

/* ---------------------------- */
/* Mumps - Constant Definitions */
/* ---------------------------- */
#ifdef GMGPOLAR_USE_MUMPS
#include "dmumps_c.h"
    /* Mumps inline functions s.t. indices match documentation */
    inline int& ICNTL(DMUMPS_STRUC_C& mumps_solver, int I) {
        return mumps_solver.icntl[(I) - 1];
    }
    inline double& CNTL(DMUMPS_STRUC_C& mumps_solver, int I) {
        return mumps_solver.cntl[(I) - 1];
    }
    inline int& INFOG(DMUMPS_STRUC_C& mumps_solver, int I) {
        return mumps_solver.infog[(I) - 1];
    }

    constexpr int USE_COMM_WORLD = -987654;
    constexpr int PAR_NOT_PARALLEL = 0;
    constexpr int PAR_PARALLEL = 1;

    constexpr int JOB_INIT = -1;
    constexpr int JOB_END = -2;
    constexpr int JOB_REMOVE_SAVED_DATA = -3;
    constexpr int JOB_FREE_INTERNAL_DATA = -4;
    constexpr int JOB_SUPPRESS_OOC_FILES = -200;

    constexpr int JOB_ANALYSIS_PHASE = 1;
    constexpr int JOB_FACTORIZATION_PHASE = 2;
    constexpr int JOB_COMPUTE_SOLUTION = 3;
    constexpr int JOB_ANALYSIS_AND_FACTORIZATION = 4;
    constexpr int JOB_FACTORIZATION_AND_SOLUTION = 5;
    constexpr int JOB_ANALYSIS_FACTORIZATION_SOLUTION = 6;
    constexpr int JOB_SAVE_INTERNAL_DATA = 7;
    constexpr int JOB_RESTORE_INTERNAL_DATA = 8;
    constexpr int JOB_DISTRIBUTE_RHS = 9;

    constexpr int SYM_UNSYMMETRIC = 0;
    constexpr int SYM_POSITIVE_DEFINITE = 1;
    constexpr int SYM_GENERAL_SYMMETRIC = 2;
#endif

// --------------------------------------- //
// Function-like macros for LIKWID markers //
// --------------------------------------- //

#ifdef GMGPOLAR_USE_LIKWID
    #include <likwid.h>
    #include <likwid-marker.h>

    #define LIKWID_INIT() LIKWID_MARKER_INIT

    #define LIKWID_REGISTER(marker)                                                                                    \
        _Pragma("omp parallel")                                                                                        \
        {                                                                                                              \
            LIKWID_MARKER_REGISTER(marker);                                                                            \
        }

    #define LIKWID_START(marker)                                                                                       \
        _Pragma("omp parallel")                                                                                        \
        {                                                                                                              \
            LIKWID_MARKER_START(marker);                                                                               \
        }

    #define LIKWID_STOP(marker)                                                                                        \
        _Pragma("omp parallel")                                                                                        \
        {                                                                                                              \
            LIKWID_MARKER_STOP(marker);                                                                                \
        }

    #define LIKWID_CLOSE() LIKWID_MARKER_CLOSE
#else
    #define LIKWID_INIT()
    #define LIKWID_REGISTER(marker)
    #define LIKWID_START(marker)
    #define LIKWID_STOP(marker)
    #define LIKWID_CLOSE()
#endif
