#pragma once

/**
 * @class DomainGeometry
 * @brief An abstract base class representing the geometric properties of a domain.
 *
 * This class provides an interface for calculating geometric transformations and their derivatives
 * for a domain in polar coordinates (r, θ). It includes pure virtual functions to compute the
 * Cartesian coordinates (Fx, Fy) from polar coordinates, as well as their partial derivatives with
 * respect to r and θ.
 *
 * Subclasses should implement the specific transformations and derivative calculations for their
 * respective geometric domains.
 */

class DomainGeometry
{
public:
    DomainGeometry()          = default;
    virtual ~DomainGeometry() = default;

    // In earlier versions denoted by 'x' and 'y'.
    virtual double Fx(const double& r, const double& theta) const = 0;
    virtual double Fy(const double& r, const double& theta) const = 0;

    // In earlier versions denoted by 'Jrr', 'Jtr', 'Jrt' and 'Jtt'.
    virtual double dFx_dr(const double& r, const double& theta) const = 0;
    virtual double dFy_dr(const double& r, const double& theta) const = 0;
    virtual double dFx_dt(const double& r, const double& theta) const = 0;
    virtual double dFy_dt(const double& r, const double& theta) const = 0;
};
