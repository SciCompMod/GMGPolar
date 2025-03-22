
#include "../InputFunctions/domainGeometry.h"

#include <cmath>

inline void compute_jacobian_elements(const DomainGeometry& domain_geometry, const double& r, const double& theta,
                                      const double& sin_theta, const double& cos_theta, const double& coeff_alpha,
                                      double& arr, double& att, double& art, double& detDF)
{
    /* Calculate the elements of the Jacobian matrix for the transformation mapping */
    /* The Jacobian matrix is: */
    /* [Jrr, Jrt] */
    /* [Jtr, Jtt] */
    const double Jrr = domain_geometry.dFx_dr(r, theta, sin_theta, cos_theta);
    const double Jtr = domain_geometry.dFy_dr(r, theta, sin_theta, cos_theta);
    const double Jrt = domain_geometry.dFx_dt(r, theta, sin_theta, cos_theta);
    const double Jtt = domain_geometry.dFy_dt(r, theta, sin_theta, cos_theta);
    /* Compute the determinant of the Jacobian matrix */
    detDF = Jrr * Jtt - Jrt * Jtr;
    /* Compute the elements of the symmetric matrix: */
    /* 0.5 * alpha * DF^{-1} * DF^{-T} * |det(DF)| */
    /* which is represented by: */
    /*  [arr, 0.5*art]  */
    /*  [0.5*atr, att]  */
    arr = 0.5 * (Jtt * Jtt + Jrt * Jrt) * coeff_alpha / fabs(detDF);
    att = 0.5 * (Jtr * Jtr + Jrr * Jrr) * coeff_alpha / fabs(detDF);
    art = (-Jtt * Jtr - Jrt * Jrr) * coeff_alpha / fabs(detDF);
    /* Note that the inverse Jacobian matrix DF^{-1} is: */
    /*  1.0 / det(DF) *  */
    /*  [Jtt, -Jrt]      */
    /*  [-Jtr, Jrr]      */
}
