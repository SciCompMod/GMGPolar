#include "../../../include/DirectSolver/DirectSolverGive/directSolverGive.h"

// clang-format off
#define COMPUTE_JACOBIAN_ELEMENTS(domain_geometry, r, theta, sin_theta, cos_theta, coeff_alpha, \
    arr, att, art, detDF) \
do { \
    /* Calculate the elements of the Jacobian matrix for the transformation mapping */ \
    /* The Jacobian matrix is: */ \
    /* [Jrr, Jrt] */ \
    /* [Jtr, Jtt] */ \
    const double Jrr = domain_geometry.dFx_dr(r, theta, sin_theta, cos_theta); \
    const double Jtr = domain_geometry.dFy_dr(r, theta, sin_theta, cos_theta); \
    const double Jrt = domain_geometry.dFx_dt(r, theta, sin_theta, cos_theta); \
    const double Jtt = domain_geometry.dFy_dt(r, theta, sin_theta, cos_theta); \
    /* Compute the determinant of the Jacobian matrix */ \
    detDF = Jrr * Jtt - Jrt * Jtr; \
    /* Compute the elements of the symmetric matrix: */ \
    /* 0.5 * alpha * DF^{-1} * DF^{-T} * |det(DF)| */ \
    /* which is represented by: */ \
    /* [arr, 0.5*art] */ \
    /* [0.5*atr, att] */ \
    arr = 0.5 * (Jtt * Jtt + Jrt * Jrt) * coeff_alpha / fabs(detDF); \
    att = 0.5 * (Jtr * Jtr + Jrr * Jrr) * coeff_alpha / fabs(detDF); \
    art = (- Jtt * Jtr - Jrt * Jrr) * coeff_alpha / fabs(detDF); \
    /* Note that the inverse Jacobian matrix DF^{-1} is: */ \
    /* 1.0 / det(DF) *   */ \
    /* [Jtt, -Jrt] */ \
    /* [-Jtr, Jrr] */ \
} while(0)
// clang-format on

/* ----------------------- */
/* Boundary Symmetry Shift */
/* ----------------------- */

void DirectSolverGive::applySymmetryShiftInnerBoundary(Vector<double>& x) const
{
    assert(DirBC_Interior_);

    const auto& sin_theta_cache = level_cache_.sin_theta();
    const auto& cos_theta_cache = level_cache_.cos_theta();

    int i_r;
    double r;
    double h1, h2, k1, k2;
    double coeff1, coeff2;
    double coeff_alpha, coeff_beta;
    double arr, att, art, detDF;

    for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++)
    {
        const double theta = grid_.theta(i_theta);
        const double sin_theta = sin_theta_cache[i_theta];
        const double cos_theta = cos_theta_cache[i_theta];
        /* -------------------------- */
        /* Node on the inner boundary */
        /* -------------------------- */
        i_r = 0;
        r = grid_.radius(i_r);

        if (level_cache_.cacheDensityProfileCoefficients())
        {
            coeff_beta = level_cache_.coeff_beta()[i_r];
        }
        else
        {
            coeff_beta = density_profile_coefficients_.beta(r);
        }
        if (!level_cache_.cacheDomainGeometry())
        {
            if (level_cache_.cacheDensityProfileCoefficients())
            {
                coeff_alpha = level_cache_.coeff_alpha()[i_r];
            }
            else
            {
                coeff_alpha = density_profile_coefficients_.alpha(r);
            }
        }

        /* Compute arr, att, art, detDF value at the current node */
        if (level_cache_.cacheDomainGeometry())
        {
            const int index = grid_.index(i_r, i_theta);
            arr = level_cache_.arr()[index];
            att = level_cache_.att()[index];
            art = level_cache_.art()[index];
            detDF = level_cache_.detDF()[index];
        }
        else
        {
            COMPUTE_JACOBIAN_ELEMENTS(domain_geometry_, r, theta, sin_theta, cos_theta, coeff_alpha, arr, att, art, detDF);
        }

        h2 = grid_.radialSpacing(i_r);
        k1 = grid_.angularSpacing(i_theta - 1);
        k2 = grid_.angularSpacing(i_theta);

        coeff2 = 0.5 * (k1 + k2) / h2;

        /* Fill x(i+1,j) */
        x[grid_.index(i_r + 1, i_theta)] -= -coeff2 * arr * x[grid_.index(i_r, i_theta)]     /* Left */
                                            + 0.25 * art * x[grid_.index(i_r, i_theta + 1)]  /* Top Left */
                                            - 0.25 * art * x[grid_.index(i_r, i_theta - 1)]; /* Bottom Left */

        /* --------------------------- */
        /* Node next to inner boundary */
        /* --------------------------- */
        i_r = 1;
        r = grid_.radius(i_r);

        if (level_cache_.cacheDensityProfileCoefficients())
        {
            coeff_beta = level_cache_.coeff_beta()[i_r];
        }
        else
        {
            coeff_beta = density_profile_coefficients_.beta(r);
        }
        if (!level_cache_.cacheDomainGeometry())
        {
            if (level_cache_.cacheDensityProfileCoefficients())
            {
                coeff_alpha = level_cache_.coeff_alpha()[i_r];
            }
            else
            {
                coeff_alpha = density_profile_coefficients_.alpha(r);
            }
        }

        /* Compute arr, att, art, detDF value at the current node */
        if (level_cache_.cacheDomainGeometry())
        {
            const int index = grid_.index(i_r, i_theta);
            arr = level_cache_.arr()[index];
            att = level_cache_.att()[index];
            art = level_cache_.art()[index];
            detDF = level_cache_.detDF()[index];
        }
        else
        {
            COMPUTE_JACOBIAN_ELEMENTS(domain_geometry_, r, theta, sin_theta, cos_theta, coeff_alpha, arr, att, art, detDF);
        }

        h1 = grid_.radialSpacing(i_r - 1);
        k1 = grid_.angularSpacing(i_theta - 1);
        k2 = grid_.angularSpacing(i_theta);

        coeff1 = 0.5 * (k1 + k2) / h1;

        /* Fill x(i,j) */
        x[grid_.index(i_r, i_theta)] -= -coeff1 * arr * x[grid_.index(i_r - 1, i_theta)]; /* Left */
        /* Fill x(i,j-1) */
        x[grid_.index(i_r, i_theta - 1)] -= +0.25 * art * x[grid_.index(i_r - 1, i_theta)]; /* Top Left */
        /* Fill x(i,j+1) */
        x[grid_.index(i_r, i_theta + 1)] -= -0.25 * art * x[grid_.index(i_r - 1, i_theta)]; /* Bottom Left */
    }
}

void DirectSolverGive::applySymmetryShiftOuterBoundary(Vector<double>& x) const
{
    const auto& sin_theta_cache = level_cache_.sin_theta();
    const auto& cos_theta_cache = level_cache_.cos_theta();

    int i_r;
    double r;
    double h1, h2, k1, k2;
    double coeff1, coeff2;
    double coeff_alpha, coeff_beta;
    double arr, att, art, detDF;

    for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++)
    {
        const double theta = grid_.theta(i_theta);
        const double sin_theta = sin_theta_cache[i_theta];
        const double cos_theta = cos_theta_cache[i_theta];
        /* --------------------------- */
        /* Node next to outer boundary */
        /* --------------------------- */
        i_r = grid_.nr() - 2;
        r = grid_.radius(i_r);

        if (level_cache_.cacheDensityProfileCoefficients())
        {
            coeff_beta = level_cache_.coeff_beta()[i_r];
        }
        else
        {
            coeff_beta = density_profile_coefficients_.beta(r);
        }
        if (!level_cache_.cacheDomainGeometry())
        {
            if (level_cache_.cacheDensityProfileCoefficients())
            {
                coeff_alpha = level_cache_.coeff_alpha()[i_r];
            }
            else
            {
                coeff_alpha = density_profile_coefficients_.alpha(r);
            }
        }

        /* Compute arr, att, art, detDF value at the current node */
        if (level_cache_.cacheDomainGeometry())
        {
            const int index = grid_.index(i_r, i_theta);
            arr = level_cache_.arr()[index];
            att = level_cache_.att()[index];
            art = level_cache_.art()[index];
            detDF = level_cache_.detDF()[index];
        }
        else
        {
            COMPUTE_JACOBIAN_ELEMENTS(domain_geometry_, r, theta, sin_theta, cos_theta, coeff_alpha, arr, att, art, detDF);
        }

        h2 = grid_.radialSpacing(i_r);
        k1 = grid_.angularSpacing(i_theta - 1);
        k2 = grid_.angularSpacing(i_theta);

        coeff2 = 0.5 * (k1 + k2) / h2;

        /* Fill result(i,j) */
        x[grid_.index(i_r, i_theta)] -= -coeff2 * arr * x[grid_.index(i_r + 1, i_theta)]; /* Right */
        /* Fill result(i,j-1) */
        x[grid_.index(i_r, i_theta - 1)] -= -0.25 * art * x[grid_.index(i_r + 1, i_theta)]; /* Top Right */
        /* Fill result(i,j+1) */
        x[grid_.index(i_r, i_theta + 1)] -= +0.25 * art * x[grid_.index(i_r + 1, i_theta)]; /* Bottom Right */

        /* -------------------------- */
        /* Node on the outer boundary */
        /* -------------------------- */
        i_r = grid_.nr() - 1;
        r = grid_.radius(i_r);

        if (level_cache_.cacheDensityProfileCoefficients())
        {
            coeff_beta = level_cache_.coeff_beta()[i_r];
        }
        else
        {
            coeff_beta = density_profile_coefficients_.beta(r);
        }
        if (!level_cache_.cacheDomainGeometry())
        {
            if (level_cache_.cacheDensityProfileCoefficients())
            {
                coeff_alpha = level_cache_.coeff_alpha()[i_r];
            }
            else
            {
                coeff_alpha = density_profile_coefficients_.alpha(r);
            }
        }

        /* Compute arr, att, art, detDF value at the current node */
        if (level_cache_.cacheDomainGeometry())
        {
            const int index = grid_.index(i_r, i_theta);
            arr = level_cache_.arr()[index];
            att = level_cache_.att()[index];
            art = level_cache_.art()[index];
            detDF = level_cache_.detDF()[index];
        }
        else
        {
            COMPUTE_JACOBIAN_ELEMENTS(domain_geometry_, r, theta, sin_theta, cos_theta, coeff_alpha, arr, att, art, detDF);
        }

        h1 = grid_.radialSpacing(i_r - 1);
        k1 = grid_.angularSpacing(i_theta - 1);
        k2 = grid_.angularSpacing(i_theta);

        coeff1 = 0.5 * (k1 + k2) / h1;

        /* Fill result(i-1,j) */
        x[grid_.index(i_r - 1, i_theta)] -= -coeff1 * arr * x[grid_.index(i_r, i_theta)]     /* Right */
                                            - 0.25 * art * x[grid_.index(i_r, i_theta + 1)]  /* Top Right */
                                            + 0.25 * art * x[grid_.index(i_r, i_theta - 1)]; /* Bottom Right */
    }
}

// clang-format off
void DirectSolverGive::applySymmetryShift(Vector<double>& x) const
{
    assert(x.size() == grid_.numberOfNodes());
    assert(grid_.nr() >= 4);

    omp_set_num_threads(num_omp_threads_);

    if (num_omp_threads_ == 1)
    {
        /* Single-threaded execution */
        if (DirBC_Interior_)
        {
            applySymmetryShiftInnerBoundary(x);
        }
        applySymmetryShiftOuterBoundary(x);
    }
    else
    {
        #pragma omp parallel sections
        {
            #pragma omp section
            {
                if (DirBC_Interior_)
                {
                    applySymmetryShiftInnerBoundary(x);
                }
            }

            #pragma omp section
            {
                applySymmetryShiftOuterBoundary(x);
            }
        }
    }
}
// clang-format on
