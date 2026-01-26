#pragma once

#include <cmath>

#include "../sourceTerm.h"

#include "../../PolarGrid/polargrid.h"

class CartesianR2_Zoni_ShafranovGeometry : public SourceTerm
{
public:
    explicit CartesianR2_Zoni_ShafranovGeometry(PolarGrid const& grid, double Rmax, double elongation_kappa,
                                                double shift_delta);
    virtual ~CartesianR2_Zoni_ShafranovGeometry() = default;

    double operator()(std::size_t i_r, std::size_t i_theta) const override;

private:
    PolarGrid const& grid_;
    const double Rmax             = 1.3;
    const double elongation_kappa = 0.3;
    const double shift_delta      = 0.2;
};
