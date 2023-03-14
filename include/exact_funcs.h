/* 
* Copyright (C) 2019-2023 The GMGPolar Development Team
*
* Authors: Philippe Leleux, Christina Schwarz, Martin J. Kühn, Carola Kruse, Ulrich Rüde
*
* Contact: 
*    Carola Kruse <Carola.Kruse@CERFACS.fr>
*    Martin J. Kuehn <Martin.Kuehn@DLR.de>
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/
#pragma once

class ExactFuncs
{
public:
    virtual ~ExactFuncs()
    {
    }
    virtual double x(double r, double theta, double kappa_eps, double delta_e, double Rmax) const                  = 0;
    virtual void x(std::vector<double> const& r, double theta, double kappa_eps, double delta_e, double Rmax,
                   std::vector<double>& sol) const                                                                 = 0;
    virtual void x(double r, std::vector<double> const& theta, double kappa_eps, double delta_e, double Rmax,
                   std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const = 0;
    virtual double y(double r, double theta, double kappa_eps, double delta_e, double Rmax) const                  = 0;
    virtual void y(std::vector<double> const& r, double theta, double kappa_eps, double delta_e, double Rmax,
                   std::vector<double>& sol) const                                                                 = 0;
    virtual void y(double r, std::vector<double> const& theta, double kappa_eps, double delta_e, double Rmax,
                   std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const = 0;
    virtual double J_rr(double r, double theta, double kappa_eps, double delta_e, double Rmax) const               = 0;
    virtual void J_rr(std::vector<double> const& r, double theta, double kappa_eps, double delta_e, double Rmax,
                      std::vector<double>& sol) const                                                              = 0;
    virtual void J_rr(double r, std::vector<double> const& theta, double kappa_eps, double delta_e, double Rmax,
                      std::vector<double>& sol, std::vector<double>& sin_theta,
                      std::vector<double>& cos_theta) const                                                        = 0;
    virtual double J_rt(double r, double theta, double kappa_eps, double delta_e, double Rmax) const               = 0;
    virtual void J_rt(std::vector<double> const& r, double theta, double kappa_eps, double delta_e, double Rmax,
                      std::vector<double>& sol) const                                                              = 0;
    virtual void J_rt(double r, std::vector<double> const& theta, double kappa_eps, double delta_e, double Rmax,
                      std::vector<double>& sol, std::vector<double>& sin_theta,
                      std::vector<double>& cos_theta) const                                                        = 0;
    virtual double J_tr(double r, double theta, double kappa_eps, double delta_e, double Rmax) const               = 0;
    virtual void J_tr(std::vector<double> const& r, double theta, double kappa_eps, double delta_e, double Rmax,
                      std::vector<double>& sol) const                                                              = 0;
    virtual void J_tr(double r, std::vector<double> const& theta, double kappa_eps, double delta_e, double Rmax,
                      std::vector<double>& sol, std::vector<double>& sin_theta,
                      std::vector<double>& cos_theta) const                                                        = 0;
    virtual double J_tt(double r, double theta, double kappa_eps, double delta_e, double Rmax) const               = 0;
    virtual void J_tt(std::vector<double> const& r, double theta, double kappa_eps, double delta_e, double Rmax,
                      std::vector<double>& sol) const                                                              = 0;
    virtual void J_tt(double r, std::vector<double> const& theta, double kappa_eps, double delta_e, double Rmax,
                      std::vector<double>& sol, std::vector<double>& sin_theta,
                      std::vector<double>& cos_theta) const                                                        = 0;
    virtual double rho_glob(double r, double theta, double kappa_eps, double delta_e, double Rmax) const           = 0;
    virtual void rho_glob(std::vector<double> const& r, double theta, double kappa_eps, double delta_e, double Rmax,
                          std::vector<double>& sol) const                                                          = 0;
    virtual void rho_glob(double r, std::vector<double> const& theta, double kappa_eps, double delta_e, double Rmax,
                          std::vector<double>& sol, std::vector<double>& sin_theta,
                          std::vector<double>& cos_theta) const                                                    = 0;
    virtual double rho_pole(double r, double theta, double kappa_eps, double delta_e, double Rmax) const           = 0;
    virtual void rho_pole(std::vector<double> const& r, double theta, double kappa_eps, double delta_e, double Rmax,
                          std::vector<double>& sol) const                                                          = 0;
    virtual void rho_pole(double r, std::vector<double> const& theta, double kappa_eps, double delta_e, double Rmax,
                          std::vector<double>& sol, std::vector<double>& sin_theta,
                          std::vector<double>& cos_theta) const                                                    = 0;
    virtual double coeffs1(double r, double Rmax) const                                                            = 0;
    virtual void coeffs1(std::vector<double> const& r, double Rmax, std::vector<double>& sol) const                = 0;
    virtual double coeffs2(double r, double Rmax) const                                                            = 0;
    virtual void coeffs2(std::vector<double> const& r, double Rmax, std::vector<double>& sol) const                = 0;
    virtual double phi_exact(double r, double theta, double kappa_eps, double delta_e, double Rmax) const          = 0;
    virtual void phi_exact(std::vector<double> const& r, double theta, double kappa_eps, double delta_e, double Rmax,
                           std::vector<double>& sol) const                                                         = 0;
    virtual void phi_exact(double r, std::vector<double> const& theta, double kappa_eps, double delta_e, double Rmax,
                           std::vector<double>& sol, std::vector<double>& sin_theta,
                           std::vector<double>& cos_theta) const                                                   = 0;
};
