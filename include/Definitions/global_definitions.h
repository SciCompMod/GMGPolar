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
