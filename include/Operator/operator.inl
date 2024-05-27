inline void Operator::arr_att_art(
    const double& r, const double& theta, const double& sin_theta, const double& cos_theta, const double& coeff_alpha, 
    double& arr, double& att, double& art, double& detDF) const
{
    // Calculate the elements of the Jacobian matrix for the transformation mapping
    // The Jacobian matrix is:
    // [Jrr, Jrt]
    // [Jtr, Jtt]
    const double Jrr = (*dFx_dr_)(r, theta, sin_theta, cos_theta);
    const double Jtr = (*dFy_dr_)(r, theta, sin_theta, cos_theta);
    const double Jrt = (*dFx_dt_)(r, theta, sin_theta, cos_theta);
    const double Jtt = (*dFy_dt_)(r, theta, sin_theta, cos_theta);
    // Compute the determinant of the Jacobian matrix
    detDF = Jrr * Jtt - Jrt * Jtr;
    // Compute the elements of the symmetric matrix:
    // 0.5 * alpha * DF^{-1} * DF^{-T} * |det(DF)|
    // which is represented as:
    // [arr, 0.5*art]
    // [0.5*atr, att]
    arr = 0.5 * (Jtt * Jtt + Jrt * Jrt) * coeff_alpha / fabs(detDF);
    att = 0.5 * (Jtr * Jtr + Jrr * Jrr) * coeff_alpha / fabs(detDF);
    art = (- Jtt * Jtr - Jrt * Jrr) * coeff_alpha / fabs(detDF);
    // Note that the inverse Jacobian matrix is:
    // 1.0 / det(DF) * 
    // [Jtt, -Jrt]
    // [-Jtr, Jrr]
}