#pragma once

enum alpha_coeff
{
    SONNENDRUCKER = 0,
    ZONI          = 1,
    ZONI_SHIFTED  = 2,
    POISSON       = 3
};

enum beta_coeff
{
    ZERO = 0,
    ALPHA_INVERSE = 1
};

enum geometry_type
{
    CIRCULAR   = 0,
    SHAFRANOV  = 1,
    TRIANGULAR = 2,
    CULHAM     = 3 
};

enum problem_type
{
    CARTESIAN_R2   = 0, 
    POLAR_R6       = 1,
    CARTESIAN_R6   = 2,
    REFINED_RADIUS = 3
};