#include "../../include/GMGPolar/gmgpolar.h"

void GMGPolar::selectTestCase(){
    /* --------------- */
    /* Domain Geometry */
    switch (geometry_) {
        case GeometryType::CIRCULAR:
            domain_geometry_ = std::make_unique<CircularGeometry>(Rmax_);
            break;

        case GeometryType::SHAFRANOV:
            domain_geometry_ = std::make_unique<ShafranovGeometry>(Rmax_, kappa_eps_, delta_e_);
            break;

        case GeometryType::CZARNY:
            domain_geometry_ = std::make_unique<CzarnyGeometry>(Rmax_, kappa_eps_, delta_e_);
            break;

        case GeometryType::CULHAM:
            domain_geometry_ = std::make_unique<CulhamGeometry>(Rmax_);
            break;

        default:
            throw std::runtime_error("Invalid geometry.\n");
    }

    /* ---------------------------- */
    /* Density Profile Coefficients */
    switch (alpha_) {
        case AlphaCoeff::POISSON:
            density_profile_coefficients_ = std::make_unique<PoissonCoefficients>(Rmax_, alpha_jump_);
            break;

        case AlphaCoeff::SONNENDRUCKER:
            switch (beta_) {
                case BetaCoeff::ZERO:
                    density_profile_coefficients_ = std::make_unique<SonnendruckerCoefficients>(Rmax_, alpha_jump_);
                    break;
                case BetaCoeff::ALPHA_INVERSE:
                    density_profile_coefficients_ = std::make_unique<SonnendruckerGyroCoefficients>(Rmax_, alpha_jump_);
                    break;
                default:
                    throw std::runtime_error("Invalid beta.\n");
            }
            break;

        case AlphaCoeff::ZONI:
            switch (beta_) {
                case BetaCoeff::ZERO:
                    density_profile_coefficients_ = std::make_unique<ZoniCoefficients>(Rmax_, alpha_jump_);
                    break;
                case BetaCoeff::ALPHA_INVERSE:
                    density_profile_coefficients_ = std::make_unique<ZoniGyroCoefficients>(Rmax_, alpha_jump_);
                    break;
                default:
                    throw std::runtime_error("Invalid beta.\n");
            }
            break;

        case AlphaCoeff::ZONI_SHIFTED:
            switch (beta_) {
                case BetaCoeff::ZERO:
                    density_profile_coefficients_ = std::make_unique<ZoniShiftedCoefficients>(Rmax_, alpha_jump_);
                    break;
                case BetaCoeff::ALPHA_INVERSE:
                    density_profile_coefficients_ = std::make_unique<ZoniShiftedGyroCoefficients>(Rmax_, alpha_jump_);
                    break;
                default:
                    throw std::runtime_error("Invalid beta.\n");
            }
            break;

        default:
            throw std::runtime_error("Invalid alpha.\n");
    }

    /* ------------------------------------ */
    /* Exact Solution & Boundary Conditions */
    switch (problem_) {
        case ProblemType::CARTESIAN_R2:
            switch (geometry_) {
                case GeometryType::CIRCULAR:
                    exact_solution_ = std::make_unique<CartesianR2_CircularGeometry>(Rmax_);
                    boundary_conditions_ = std::make_unique<CartesianR2_Boundary_CircularGeometry>(Rmax_);
                    break;
                case GeometryType::SHAFRANOV:
                    exact_solution_ = std::make_unique<CartesianR2_ShafranovGeometry>(Rmax_, kappa_eps_, delta_e_);
                    boundary_conditions_ = std::make_unique<CartesianR2_Boundary_ShafranovGeometry>(Rmax_, kappa_eps_, delta_e_);
                    break;
                case GeometryType::CZARNY:
                    exact_solution_ = std::make_unique<CartesianR2_CzarnyGeometry>(Rmax_, kappa_eps_, delta_e_);
                    boundary_conditions_ = std::make_unique<CartesianR2_Boundary_CzarnyGeometry>(Rmax_, kappa_eps_, delta_e_);
                    break;
                default:
                    throw std::runtime_error("Invalid geometry for configuration.\n");
            }
            break;

        case ProblemType::CARTESIAN_R6:
            switch (geometry_) {
                case GeometryType::CIRCULAR:
                    exact_solution_ = std::make_unique<CartesianR6_CircularGeometry>(Rmax_);
                    boundary_conditions_ = std::make_unique<CartesianR6_Boundary_CircularGeometry>(Rmax_);
                    break;
                case GeometryType::SHAFRANOV:
                    exact_solution_ = std::make_unique<CartesianR6_ShafranovGeometry>(Rmax_, kappa_eps_, delta_e_);
                    boundary_conditions_ = std::make_unique<CartesianR6_Boundary_ShafranovGeometry>(Rmax_, kappa_eps_, delta_e_);
                    break;
                case GeometryType::CZARNY:
                    exact_solution_ = std::make_unique<CartesianR6_CzarnyGeometry>(Rmax_, kappa_eps_, delta_e_);
                    boundary_conditions_ = std::make_unique<CartesianR6_Boundary_CzarnyGeometry>(Rmax_, kappa_eps_, delta_e_);
                    break;
                default:
                    throw std::runtime_error("Invalid geometry for configuration.\n");
            }
            break;

        case ProblemType::POLAR_R6:
            switch (geometry_) {
                case GeometryType::CIRCULAR:
                    exact_solution_ = std::make_unique<PolarR6_CircularGeometry>(Rmax_);
                    boundary_conditions_ = std::make_unique<PolarR6_Boundary_CircularGeometry>(Rmax_);
                    break;
                case GeometryType::SHAFRANOV:
                    exact_solution_ = std::make_unique<PolarR6_ShafranovGeometry>(Rmax_, kappa_eps_, delta_e_);
                    boundary_conditions_ = std::make_unique<PolarR6_Boundary_ShafranovGeometry>(Rmax_, kappa_eps_, delta_e_);
                    break;
                case GeometryType::CZARNY:
                    exact_solution_ = std::make_unique<PolarR6_CzarnyGeometry>(Rmax_, kappa_eps_, delta_e_);
                    boundary_conditions_ = std::make_unique<PolarR6_Boundary_CzarnyGeometry>(Rmax_, kappa_eps_, delta_e_);
                    break;
                case GeometryType::CULHAM:
                    exact_solution_ = std::make_unique<PolarR6_CulhamGeometry>(Rmax_);
                    boundary_conditions_ = std::make_unique<PolarR6_Boundary_CulhamGeometry>(Rmax_);
                    break;
                default:
                    throw std::runtime_error("Invalid geometry for configuration.\n");
            }
            break;
            
        case ProblemType::REFINED_RADIUS:
            switch (geometry_) {
                case GeometryType::CIRCULAR:
                    exact_solution_ = std::make_unique<Refined_CircularGeometry>(Rmax_);
                    boundary_conditions_ = std::make_unique<Refined_Boundary_CircularGeometry>(Rmax_);
                    break;
                case GeometryType::SHAFRANOV:
                    exact_solution_ = std::make_unique<Refined_ShafranovGeometry>(Rmax_, kappa_eps_, delta_e_);
                    boundary_conditions_ = std::make_unique<Refined_Boundary_ShafranovGeometry>(Rmax_, kappa_eps_, delta_e_);
                    break;
                case GeometryType::CZARNY:
                    exact_solution_ = std::make_unique<Refined_CzarnyGeometry>(Rmax_, kappa_eps_, delta_e_);
                    boundary_conditions_ = std::make_unique<Refined_Boundary_CzarnyGeometry>(Rmax_, kappa_eps_, delta_e_);
                    break;
                case GeometryType::CULHAM:
                    exact_solution_ = std::make_unique<Refined_CulhamGeometry>(Rmax_);
                    boundary_conditions_ = std::make_unique<Refined_Boundary_CulhamGeometry>(Rmax_);
                    break;
                default:
                    throw std::runtime_error("Invalid geometry for configuration.\n");
            }
            break;

        default:
            throw std::runtime_error("Invalid problem.\n");
    }

    /* ------------------- */
    /* Source Term (rhs_f) */
    switch (problem_) {
        case ProblemType::CARTESIAN_R2:

            switch (geometry_) {
                case GeometryType::CIRCULAR:

                    switch (alpha_) {
                        case AlphaCoeff::POISSON:
                            source_term_ = std::make_unique<CartesianR2_Poisson_CircularGeometry>(Rmax_);
                            break;
                        case AlphaCoeff::SONNENDRUCKER:
                            switch (beta_) {
                                case BetaCoeff::ZERO:
                                    source_term_ = std::make_unique<CartesianR2_Sonnendrucker_CircularGeometry>(Rmax_);
                                    break;
                                case BetaCoeff::ALPHA_INVERSE:
                                    source_term_ = std::make_unique<CartesianR2_SonnendruckerGyro_CircularGeometry>(Rmax_);
                                    break;
                                default:
                                    throw std::runtime_error("Invalid beta.\n");
                                
                            }
                            break;
                        case AlphaCoeff::ZONI:
                            switch (beta_) {
                                case BetaCoeff::ZERO:
                                    source_term_ = std::make_unique<CartesianR2_Zoni_CircularGeometry>(Rmax_);
                                    break;
                                case BetaCoeff::ALPHA_INVERSE:
                                    source_term_ = std::make_unique<CartesianR2_ZoniGyro_CircularGeometry>(Rmax_);
                                    break;
                                default:
                                    throw std::runtime_error("Invalid beta.\n");
                            }
                            break;
                        case AlphaCoeff::ZONI_SHIFTED:
                            switch (beta_) {
                                case BetaCoeff::ZERO:
                                    source_term_ = std::make_unique<CartesianR2_ZoniShifted_CircularGeometry>(Rmax_);
                                    break;
                                case BetaCoeff::ALPHA_INVERSE:
                                    source_term_ = std::make_unique<CartesianR2_ZoniShiftedGyro_CircularGeometry>(Rmax_);
                                    break;
                                default:
                                    throw std::runtime_error("Invalid beta.\n");
                            }
                            break;
                        default:
                            throw std::runtime_error("Invalid alpha.\n");
                    }
                    break;

                case GeometryType::SHAFRANOV:

                    switch (alpha_) {
                        case AlphaCoeff::POISSON:
                            source_term_ = std::make_unique<CartesianR2_Poisson_ShafranovGeometry>(Rmax_, kappa_eps_, delta_e_);
                            break;
                        case AlphaCoeff::SONNENDRUCKER:
                            switch (beta_) {
                                case BetaCoeff::ZERO:
                                    source_term_ = std::make_unique<CartesianR2_Sonnendrucker_ShafranovGeometry>(Rmax_, kappa_eps_, delta_e_);
                                    break;
                                case BetaCoeff::ALPHA_INVERSE:
                                    source_term_ = std::make_unique<CartesianR2_SonnendruckerGyro_ShafranovGeometry>(Rmax_, kappa_eps_, delta_e_);
                                    break;
                                default:
                                    throw std::runtime_error("Invalid beta.\n");
                            }
                            break;
                        case AlphaCoeff::ZONI:
                            switch (beta_) {
                                case BetaCoeff::ZERO:
                                    source_term_ = std::make_unique<CartesianR2_Zoni_ShafranovGeometry>(Rmax_, kappa_eps_, delta_e_);
                                    break;
                                case BetaCoeff::ALPHA_INVERSE:
                                    source_term_ = std::make_unique<CartesianR2_ZoniGyro_ShafranovGeometry>(Rmax_, kappa_eps_, delta_e_);
                                    break;
                                default:
                                    throw std::runtime_error("Invalid beta.\n");
                            }
                            break;
                        case AlphaCoeff::ZONI_SHIFTED:
                            switch (beta_) {
                                case BetaCoeff::ZERO:
                                    source_term_ = std::make_unique<CartesianR2_ZoniShifted_ShafranovGeometry>(Rmax_, kappa_eps_, delta_e_);
                                    break;
                                case BetaCoeff::ALPHA_INVERSE:
                                    source_term_ = std::make_unique<CartesianR2_ZoniShiftedGyro_ShafranovGeometry>(Rmax_, kappa_eps_, delta_e_);
                                    break;
                                default:
                                    throw std::runtime_error("Invalid beta.\n");
                            }
                            break;
                        default:
                            throw std::runtime_error("Invalid alpha.\n");
                    }
                    break;

                case GeometryType::CZARNY:

                    switch (alpha_) {
                        case AlphaCoeff::POISSON:
                            source_term_ = std::make_unique<CartesianR2_Poisson_CzarnyGeometry>(Rmax_, kappa_eps_, delta_e_);
                            break;
                        case AlphaCoeff::SONNENDRUCKER:
                            switch (beta_) {
                                case BetaCoeff::ZERO:
                                    source_term_ = std::make_unique<CartesianR2_Sonnendrucker_CzarnyGeometry>(Rmax_, kappa_eps_, delta_e_);
                                    break;
                                case BetaCoeff::ALPHA_INVERSE:
                                    source_term_ = std::make_unique<CartesianR2_SonnendruckerGyro_CzarnyGeometry>(Rmax_, kappa_eps_, delta_e_);
                                    break;
                                default:
                                    throw std::runtime_error("Invalid beta.\n");
                            }
                            break;
                        case AlphaCoeff::ZONI:
                            switch (beta_) {
                                case BetaCoeff::ZERO:
                                    source_term_ = std::make_unique<CartesianR2_Zoni_CzarnyGeometry>(Rmax_, kappa_eps_, delta_e_);
                                    break;
                                case BetaCoeff::ALPHA_INVERSE:
                                    source_term_ = std::make_unique<CartesianR2_ZoniGyro_CzarnyGeometry>(Rmax_, kappa_eps_, delta_e_);
                                    break;
                                default:
                                    throw std::runtime_error("Invalid beta.\n");
                            }
                            break;
                        case AlphaCoeff::ZONI_SHIFTED:
                            switch (beta_) {
                                case BetaCoeff::ZERO:
                                    source_term_ = std::make_unique<CartesianR2_ZoniShifted_CzarnyGeometry>(Rmax_, kappa_eps_, delta_e_);
                                    break;
                                case BetaCoeff::ALPHA_INVERSE:
                                    source_term_ = std::make_unique<CartesianR2_ZoniShiftedGyro_CzarnyGeometry>(Rmax_, kappa_eps_, delta_e_);
                                    break;
                                default:
                                    throw std::runtime_error("Invalid beta.\n");
                            }
                            break;
                        default:
                            throw std::runtime_error("Invalid alpha.\n");
                    }
                    break;

                default:
                    throw std::runtime_error("Invalid geometry for configuration.\n");
            }
            break;

        case ProblemType::CARTESIAN_R6:

            switch (geometry_) {
                case GeometryType::CIRCULAR:

                    switch (alpha_) {
                        case AlphaCoeff::POISSON:
                            source_term_ = std::make_unique<CartesianR6_Poisson_CircularGeometry>(Rmax_);
                            break;
                        case AlphaCoeff::SONNENDRUCKER:
                            switch (beta_) {
                                case BetaCoeff::ZERO:
                                    source_term_ = std::make_unique<CartesianR6_Sonnendrucker_CircularGeometry>(Rmax_);
                                    break;
                                case BetaCoeff::ALPHA_INVERSE:
                                    source_term_ = std::make_unique<CartesianR6_SonnendruckerGyro_CircularGeometry>(Rmax_);
                                    break;
                                default:
                                    throw std::runtime_error("Invalid beta.\n");
                            }
                            break;
                        case AlphaCoeff::ZONI:
                            switch (beta_) {
                                case BetaCoeff::ZERO:
                                    source_term_ = std::make_unique<CartesianR6_Zoni_CircularGeometry>(Rmax_);
                                    break;
                                case BetaCoeff::ALPHA_INVERSE:
                                    source_term_ = std::make_unique<CartesianR6_ZoniGyro_CircularGeometry>(Rmax_);
                                    break;
                                default:
                                    throw std::runtime_error("Invalid beta.\n");
                            }
                            break;
                        case AlphaCoeff::ZONI_SHIFTED:
                            switch (beta_) {
                                case BetaCoeff::ZERO:
                                    source_term_ = std::make_unique<CartesianR6_ZoniShifted_CircularGeometry>(Rmax_);
                                    break;
                                case BetaCoeff::ALPHA_INVERSE:
                                    source_term_ = std::make_unique<CartesianR6_ZoniShiftedGyro_CircularGeometry>(Rmax_);
                                    break;
                                default:
                                    throw std::runtime_error("Invalid beta.\n");
                            }
                            break;
                        default:
                            throw std::runtime_error("Invalid alpha.\n");
                    }
                    break;

                case GeometryType::SHAFRANOV:

                    switch (alpha_) {
                        case AlphaCoeff::POISSON:
                            source_term_ = std::make_unique<CartesianR6_Poisson_ShafranovGeometry>(Rmax_, kappa_eps_, delta_e_);
                            break;
                        case AlphaCoeff::SONNENDRUCKER:
                            switch (beta_) {
                                case BetaCoeff::ZERO:
                                    source_term_ = std::make_unique<CartesianR6_Sonnendrucker_ShafranovGeometry>(Rmax_, kappa_eps_, delta_e_);
                                    break;
                                case BetaCoeff::ALPHA_INVERSE:
                                    source_term_ = std::make_unique<CartesianR6_SonnendruckerGyro_ShafranovGeometry>(Rmax_, kappa_eps_, delta_e_);
                                    break;
                                default:
                                    throw std::runtime_error("Invalid beta.\n");
                            }
                            break;
                        case AlphaCoeff::ZONI:
                            switch (beta_) {
                                case BetaCoeff::ZERO:
                                    source_term_ = std::make_unique<CartesianR6_Zoni_ShafranovGeometry>(Rmax_, kappa_eps_, delta_e_);
                                    break;
                                case BetaCoeff::ALPHA_INVERSE:
                                    source_term_ = std::make_unique<CartesianR6_ZoniGyro_ShafranovGeometry>(Rmax_, kappa_eps_, delta_e_);
                                    break;
                                default:
                                    throw std::runtime_error("Invalid beta.\n");
                            }
                            break;
                        case AlphaCoeff::ZONI_SHIFTED:
                            switch (beta_) {
                                case BetaCoeff::ZERO:
                                    source_term_ = std::make_unique<CartesianR6_ZoniShifted_ShafranovGeometry>(Rmax_, kappa_eps_, delta_e_);
                                    break;
                                case BetaCoeff::ALPHA_INVERSE:
                                    source_term_ = std::make_unique<CartesianR6_ZoniShiftedGyro_ShafranovGeometry>(Rmax_, kappa_eps_, delta_e_);
                                    break;
                                default:
                                    throw std::runtime_error("Invalid beta.\n");
                            }
                            break;
                        default:
                            throw std::runtime_error("Invalid alpha.\n");
                    }
                    break;

                case GeometryType::CZARNY:

                    switch (alpha_) {
                        case AlphaCoeff::POISSON:
                            source_term_ = std::make_unique<CartesianR6_Poisson_CzarnyGeometry>(Rmax_, kappa_eps_, delta_e_);
                            break;
                        case AlphaCoeff::SONNENDRUCKER:
                            switch (beta_) {
                                case BetaCoeff::ZERO:
                                    source_term_ = std::make_unique<CartesianR6_Sonnendrucker_CzarnyGeometry>(Rmax_, kappa_eps_, delta_e_);
                                    break;
                                case BetaCoeff::ALPHA_INVERSE:
                                    source_term_ = std::make_unique<CartesianR6_SonnendruckerGyro_CzarnyGeometry>(Rmax_, kappa_eps_, delta_e_);
                                    break;
                                default:
                                    throw std::runtime_error("Invalid beta.\n");
                            }
                            break;
                        case AlphaCoeff::ZONI:
                            switch (beta_) {
                                case BetaCoeff::ZERO:
                                    source_term_ = std::make_unique<CartesianR6_Zoni_CzarnyGeometry>(Rmax_, kappa_eps_, delta_e_);
                                    break;
                                case BetaCoeff::ALPHA_INVERSE:
                                    source_term_ = std::make_unique<CartesianR6_ZoniGyro_CzarnyGeometry>(Rmax_, kappa_eps_, delta_e_);
                                    break;
                                default:
                                    throw std::runtime_error("Invalid beta.\n");
                            }
                            break;
                        case AlphaCoeff::ZONI_SHIFTED:
                            switch (beta_) {
                                case BetaCoeff::ZERO:
                                    source_term_ = std::make_unique<CartesianR6_ZoniShifted_CzarnyGeometry>(Rmax_, kappa_eps_, delta_e_);
                                    break;
                                case BetaCoeff::ALPHA_INVERSE:
                                    source_term_ = std::make_unique<CartesianR6_ZoniShiftedGyro_CzarnyGeometry>(Rmax_, kappa_eps_, delta_e_);
                                    break;
                                default:
                                    throw std::runtime_error("Invalid beta.\n");
                            }
                            break;
                        default:
                            throw std::runtime_error("Invalid alpha.\n");
                    }
                    break;

                default:
                    throw std::runtime_error("Invalid geometry for configuration.\n");
            }
            break;

        case ProblemType::POLAR_R6:

            switch (geometry_) {
                case GeometryType::CIRCULAR:

                    switch (alpha_) {
                        case AlphaCoeff::POISSON:
                            source_term_ = std::make_unique<PolarR6_Poisson_CircularGeometry>(Rmax_);
                            break;
                        case AlphaCoeff::SONNENDRUCKER:
                            switch (beta_) {
                                case BetaCoeff::ZERO:
                                    source_term_ = std::make_unique<PolarR6_Sonnendrucker_CircularGeometry>(Rmax_);
                                    break;
                                case BetaCoeff::ALPHA_INVERSE:
                                    source_term_ = std::make_unique<PolarR6_SonnendruckerGyro_CircularGeometry>(Rmax_);
                                    break;
                                default:
                                    throw std::runtime_error("Invalid beta.\n");
                            }
                            break;
                        case AlphaCoeff::ZONI:
                            switch (beta_) {
                                case BetaCoeff::ZERO:
                                    source_term_ = std::make_unique<PolarR6_Zoni_CircularGeometry>(Rmax_);
                                    break;
                                case BetaCoeff::ALPHA_INVERSE:
                                    source_term_ = std::make_unique<PolarR6_ZoniGyro_CircularGeometry>(Rmax_);
                                    break;
                                default:
                                    throw std::runtime_error("Invalid beta.\n");
                            }
                            break;
                        case AlphaCoeff::ZONI_SHIFTED:
                            switch (beta_) {
                                case BetaCoeff::ZERO:
                                    source_term_ = std::make_unique<PolarR6_ZoniShifted_CircularGeometry>(Rmax_);
                                    break;
                                case BetaCoeff::ALPHA_INVERSE:
                                    source_term_ = std::make_unique<PolarR6_ZoniShiftedGyro_CircularGeometry>(Rmax_);
                                    break;
                                default:
                                    throw std::runtime_error("Invalid beta.\n");
                            }
                            break;
                        default:
                            throw std::runtime_error("Invalid alpha.\n");
                    }
                    break;

                case GeometryType::SHAFRANOV:

                    switch (alpha_) {
                        case AlphaCoeff::POISSON:
                            source_term_ = std::make_unique<PolarR6_Poisson_ShafranovGeometry>(Rmax_, kappa_eps_, delta_e_);
                            break;
                        case AlphaCoeff::SONNENDRUCKER:
                            switch (beta_) {
                                case BetaCoeff::ZERO:
                                    source_term_ = std::make_unique<PolarR6_Sonnendrucker_ShafranovGeometry>(Rmax_, kappa_eps_, delta_e_);
                                    break;
                                case BetaCoeff::ALPHA_INVERSE:
                                    source_term_ = std::make_unique<PolarR6_SonnendruckerGyro_ShafranovGeometry>(Rmax_, kappa_eps_, delta_e_);
                                    break;
                                default:
                                    throw std::runtime_error("Invalid beta.\n");
                            }
                            break;
                        case AlphaCoeff::ZONI:
                            switch (beta_) {
                                case BetaCoeff::ZERO:
                                    source_term_ = std::make_unique<PolarR6_Zoni_ShafranovGeometry>(Rmax_, kappa_eps_, delta_e_);
                                    break;
                                case BetaCoeff::ALPHA_INVERSE:
                                    source_term_ = std::make_unique<PolarR6_ZoniGyro_ShafranovGeometry>(Rmax_, kappa_eps_, delta_e_);
                                    break;
                                default:
                                    throw std::runtime_error("Invalid beta.\n");
                            }
                            break;
                        case AlphaCoeff::ZONI_SHIFTED:
                            switch (beta_) {
                                case BetaCoeff::ZERO:
                                    source_term_ = std::make_unique<PolarR6_ZoniShifted_ShafranovGeometry>(Rmax_, kappa_eps_, delta_e_);
                                    break;
                                case BetaCoeff::ALPHA_INVERSE:
                                    source_term_ = std::make_unique<PolarR6_ZoniShiftedGyro_ShafranovGeometry>(Rmax_, kappa_eps_, delta_e_);
                                    break;
                                default:
                                    throw std::runtime_error("Invalid beta.\n");
                            }
                            break;
                        default:
                            throw std::runtime_error("Invalid alpha.\n");
                    }
                    break;

                case GeometryType::CZARNY:

                    switch (alpha_) {
                        case AlphaCoeff::POISSON:
                            source_term_ = std::make_unique<PolarR6_Poisson_CzarnyGeometry>(Rmax_, kappa_eps_, delta_e_);
                            break;
                        case AlphaCoeff::SONNENDRUCKER:
                            switch (beta_) {
                                case BetaCoeff::ZERO:
                                    source_term_ = std::make_unique<PolarR6_Sonnendrucker_CzarnyGeometry>(Rmax_, kappa_eps_, delta_e_);
                                    break;
                                case BetaCoeff::ALPHA_INVERSE:
                                    source_term_ = std::make_unique<PolarR6_SonnendruckerGyro_CzarnyGeometry>(Rmax_, kappa_eps_, delta_e_);
                                    break;
                                default:
                                    throw std::runtime_error("Invalid beta.\n");
                            }
                            break;
                        case AlphaCoeff::ZONI:
                            switch (beta_) {
                                case BetaCoeff::ZERO:
                                    source_term_ = std::make_unique<PolarR6_Zoni_CzarnyGeometry>(Rmax_, kappa_eps_, delta_e_);
                                    break;
                                case BetaCoeff::ALPHA_INVERSE:
                                    source_term_ = std::make_unique<PolarR6_ZoniGyro_CzarnyGeometry>(Rmax_, kappa_eps_, delta_e_);
                                    break;
                                default:
                                    throw std::runtime_error("Invalid beta.\n");
                            }
                            break;
                        case AlphaCoeff::ZONI_SHIFTED:
                            switch (beta_) {
                                case BetaCoeff::ZERO:
                                    source_term_ = std::make_unique<PolarR6_ZoniShifted_CzarnyGeometry>(Rmax_, kappa_eps_, delta_e_);
                                    break;
                                case BetaCoeff::ALPHA_INVERSE:
                                    source_term_ = std::make_unique<PolarR6_ZoniShiftedGyro_CzarnyGeometry>(Rmax_, kappa_eps_, delta_e_);
                                    break;
                                default:
                                    throw std::runtime_error("Invalid beta.\n");
                            }
                            break;
                        default:
                            throw std::runtime_error("Invalid alpha.\n");
                    }
                    break;

                case GeometryType::CULHAM:
                    switch (alpha_) {
                        case AlphaCoeff::ZONI_SHIFTED:
                            switch (beta_) {
                                case BetaCoeff::ALPHA_INVERSE:
                                    source_term_ = std::make_unique<PolarR6_ZoniShiftedGyro_CulhamGeometry>(Rmax_);
                                    break;
                                default:
                                    throw std::runtime_error("Invalid beta for configuration.\n");
                            }
                            break;
                        default:
                            throw std::runtime_error("Invalid alpha for configuration.\n");
                    }
                    break;

                default:
                    throw std::runtime_error("Invalid geometry for configuration.\n");
            }
            break;


        case ProblemType::REFINED_RADIUS:

            switch (geometry_) {
                case GeometryType::CIRCULAR:

                    switch (alpha_) {
                        case AlphaCoeff::ZONI_SHIFTED:
                            switch (beta_) {
                                case BetaCoeff::ALPHA_INVERSE:
                                    source_term_ = std::make_unique<Refined_ZoniShiftedGyro_CircularGeometry>(Rmax_);
                                    break;
                                default:
                                    throw std::runtime_error("Invalid beta for configuration.\n");
                            }
                            break;
                        default:
                            throw std::runtime_error("Invalid alpha for configuration.\n");
                    }
                    break;

                case GeometryType::SHAFRANOV:

                    switch (alpha_) {
                        case AlphaCoeff::ZONI_SHIFTED:
                            switch (beta_) {
                                case BetaCoeff::ALPHA_INVERSE:
                                    source_term_ = std::make_unique<Refined_ZoniShiftedGyro_ShafranovGeometry>(Rmax_, kappa_eps_, delta_e_);
                                    break;
                                default:
                                    throw std::runtime_error("Invalid beta for configuration.\n");
                            }
                            break;
                        default:
                            throw std::runtime_error("Invalid alpha for configuration.\n");
                    }
                    break;

                case GeometryType::CZARNY:

                    switch (alpha_) {
                        case AlphaCoeff::ZONI_SHIFTED:
                            switch (beta_) {
                                case BetaCoeff::ALPHA_INVERSE:
                                    source_term_ = std::make_unique<Refined_ZoniShiftedGyro_CzarnyGeometry>(Rmax_, kappa_eps_, delta_e_);
                                    break;
                                default:
                                    throw std::runtime_error("Invalid beta for configuration.\n");
                            }
                            break;
                        default:
                            throw std::runtime_error("Invalid alpha for configuration.\n"); 
                    }
                    break;

                case GeometryType::CULHAM:

                    switch (alpha_) {
                        case AlphaCoeff::ZONI_SHIFTED:
                            switch (beta_) {
                                case BetaCoeff::ALPHA_INVERSE:
                                    source_term_ = std::make_unique<Refined_ZoniShiftedGyro_CulhamGeometry>(Rmax_);
                                    break;
                                default:
                                    throw std::runtime_error("Invalid beta for configuration.\n");
                            }
                            break;
                        default:
                            throw std::runtime_error("Invalid alpha for configuration.\n"); 
                    }
                    break;

                default:
                    throw std::runtime_error("Invalid geometry for configuration.\n");
            }
            break;

        default:
            throw std::runtime_error("Invalid problem.\n");
    }
}
