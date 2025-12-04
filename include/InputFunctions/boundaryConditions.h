#pragma once

class BoundaryConditions
{
public:
    BoundaryConditions()          = default;
    virtual ~BoundaryConditions() = default;

    virtual double u_D(double r, double theta) const = 0;
    // Only used if DirBC_Interior = true
    virtual double u_D_Interior(double r, double theta) const = 0;
};
