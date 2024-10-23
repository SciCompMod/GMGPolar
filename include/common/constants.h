#pragma once

/* ---------------------------------- */
/* GMGPolar - Enumeration Definitions */
/* ---------------------------------- */

enum class ImplementationType
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
    EUCLIDEAN = 0,             // Corresponds to the L2 norm
    WEIGHTED_EUCLIDEAN = 1,    // A weighted version of the L2 norm
    INFINITY_NORM = 2          // Corresponds to the L∞ norm
};

enum class ExtrapolationType
{
    NONE = 0,
    IMPLICIT_EXTRAPOLATION = 1,
    IMPLICIT_FULL_GRID_SMOOTHING = 2,
    COMBINED = 3,
};

/* Smoother Colors */
enum class SmootherColor {
    Black = 0,
    White = 1,
};

/* -----------*/
/* Test Cases */
/* -----------*/

/* Geometry Types - domain_geometry */
enum class GeometryType
{
    CIRCULAR = 0,
    SHAFRANOV = 1,
    CZARNY = 2,
    CULHAM = 3
};

/* Test Problem Types - exact_solution */
enum class ProblemType
{
    CARTESIAN_R2 = 0,
    CARTESIAN_R6 = 1,
    POLAR_R6 = 2,
    REFINED_RADIUS = 3
};

/* Alpha Coefficient Types - profile_coefficient alpha */
enum class AlphaCoeff
{
    POISSON = 0,
    SONNENDRUCKER = 1,
    ZONI = 2,
    ZONI_SHIFTED = 3
};

/* Beta Coefficient Types - profile_coefficient beta */
enum class BetaCoeff
{
    ZERO = 0,
    ALPHA_INVERSE = 1
};

/* ---------------------------- */
/* Mumps - Constant Definitions */
/* ---------------------------- */

/* Mumps macro s.t. indices match documentation */
#define ICNTL(I) icntl[(I)-1]
#define CNTL(I) cntl[(I)-1]

#define USE_COMM_WORLD -987654
#define PAR_NOT_PARALLEL 0
#define PAR_PARALLEL 1

#define JOB_INIT -1
#define JOB_END -2
#define JOB_REMOVE_SAVED_DATA -3
#define JOB_FREE_INTERNAL_DATA -4
#define JOB_SUPPRESS_OOC_FILES -200

#define JOB_ANALYSIS_PHASE 1
#define JOB_FACTORIZATION_PHASE 2
#define JOB_COMPUTE_SOLUTION 3
#define JOB_ANALYSIS_AND_FACTORIZATION 4
#define JOB_FACTORIZATION_AND_SOLUTION 5
#define JOB_ANALYSIS_FACTORIZATION_SOLUTION 6
#define JOB_SAVE_INTERNAL_DATA 7
#define JOB_RESTORE_INTERNAL_DATA 8
#define JOB_DISTRIBUTE_RHS 9

#define SYM_UNSYMMETRIC 0
#define SYM_POSITIVE_DEFINITE 1
#define SYM_GENERAL_SYMMETRIC 2