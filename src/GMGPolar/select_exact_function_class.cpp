#include "../../include/GMGPolar/gmgpolar.h"

std::shared_ptr<ExactFunctions> GMGPolar::selectExactFunctionsClass(){
    // return std::make_shared<CartesianR2SonnendruckerCircular>();
    switch (beta) {
    case beta_coeff::ZERO:
        switch (alpha) {
        case alpha_coeff::SONNENDRUCKER:
            switch (problem) {
            case problem_type::CARTESIAN_R2:
                switch (geometry) {
                case geometry_type::CIRCULAR:
                    std::cout << "Using CartesianR2SonnendruckerCircular\n";
                    return std::make_shared<CartesianR2SonnendruckerCircular>();
                case geometry_type::SHAFRANOV:
                    std::cout << "Using CartesianR2SonnendruckerShafranov\n";
                    return std::make_shared<CartesianR2SonnendruckerShafranov>();
                case geometry_type::TRIANGULAR:
                    std::cout << "Using CartesianR2SonnendruckerTriangular\n";
                    return std::make_shared<CartesianR2SonnendruckerTriangular>();
                default:
                    throw std::runtime_error("Wrong configuration for the geometry\n");
                }
            case problem_type::POLAR_R6:
                switch (geometry) {
                case geometry_type::CIRCULAR:
                    std::cout << "Using PolarR6SonnendruckerCircular\n";
                    return std::make_shared<PolarR6SonnendruckerCircular>();
                case geometry_type::SHAFRANOV:
                    std::cout << "Using PolarR6SonnendruckerShafranov\n";
                    return std::make_shared<PolarR6SonnendruckerShafranov>();
                case geometry_type::TRIANGULAR:
                    std::cout << "Using PolarR6SonnendruckerTriangular\n";
                    return std::make_shared<PolarR6SonnendruckerTriangular>();
                default:
                    throw std::runtime_error("Wrong configuration for the geometry\n");
                }
            case problem_type::CARTESIAN_R6:
                switch (geometry) {
                case geometry_type::CIRCULAR:
                    std::cout << "Using CartesianR6SonnendruckerCircular\n";
                    return std::make_shared<CartesianR6SonnendruckerCircular>();
                case geometry_type::SHAFRANOV:
                    std::cout << "Using CartesianR6SonnendruckerShafranov\n";
                    return std::make_shared<CartesianR6SonnendruckerShafranov>();
                case geometry_type::TRIANGULAR:
                    std::cout << "Using CartesianR6SonnendruckerTriangular\n";
                    return std::make_shared<CartesianR6SonnendruckerTriangular>();
                default:
                    throw std::runtime_error("Wrong configuration for the geometry\n");
                }
            default:
                throw std::runtime_error("Wrong configuration for the problem\n");  
            }
        case alpha_coeff::ZONI:
            switch (problem) {
            case problem_type::CARTESIAN_R2:
                switch (geometry) {
                case geometry_type::CIRCULAR:
                    std::cout << "Using CartesianR2ZoniCircular\n";
                    return std::make_shared<CartesianR2ZoniCircular>();
                case geometry_type::SHAFRANOV:
                    std::cout << "Using CartesianR2ZoniShafranov\n";
                    return std::make_shared<CartesianR2ZoniShafranov>();
                case geometry_type::TRIANGULAR:
                    std::cout << "Using CartesianR2ZoniTriangular\n";
                    return std::make_shared<CartesianR2ZoniTriangular>();
                default:
                    throw std::runtime_error("Wrong configuration for the geometry\n");
                }
            case problem_type::POLAR_R6:
                switch (geometry) {
                case geometry_type::CIRCULAR:
                    std::cout << "Using PolarR6ZoniCircular\n";
                    return std::make_shared<PolarR6ZoniCircular>();
                case geometry_type::SHAFRANOV:
                    std::cout << "Using PolarR6ZoniShafranov\n";
                    return std::make_shared<PolarR6ZoniShafranov>();
                case geometry_type::TRIANGULAR:
                    std::cout << "Using PolarR6ZoniTriangular\n";
                    return std::make_shared<PolarR6ZoniTriangular>();
                default:
                    throw std::runtime_error("Wrong configuration for the geometry\n");
                }
            case problem_type::CARTESIAN_R6:
                switch (geometry) {
                case geometry_type::CIRCULAR:
                    std::cout << "Using CartesianR6ZoniCircular\n";
                    return std::make_shared<CartesianR6ZoniCircular>();
                case geometry_type::SHAFRANOV:
                    std::cout << "Using CartesianR6ZoniShafranov\n";
                    return std::make_shared<CartesianR6ZoniShafranov>();
                case geometry_type::TRIANGULAR:
                    std::cout << "Using CartesianR6ZoniTriangular\n";
                    return std::make_shared<CartesianR6ZoniTriangular>();
                default:
                    throw std::runtime_error("Wrong configuration for the geometry\n");
                }
            default:
                throw std::runtime_error("Wrong configuration for the problem\n");
            }
        case alpha_coeff::ZONI_SHIFTED:
            switch (problem) {
            case problem_type::CARTESIAN_R2:
                switch (geometry) {
                case geometry_type::CIRCULAR:
                    std::cout << "Using CartesianR2ZoniShiftedCircular\n";
                    return std::make_shared<CartesianR2ZoniShiftedCircular>();
                case geometry_type::SHAFRANOV:
                    std::cout << "Using CartesianR2ZoniShiftedShafranov\n";
                    return std::make_shared<CartesianR2ZoniShiftedShafranov>();
                case geometry_type::TRIANGULAR:
                    std::cout << "Using CartesianR2ZoniShiftedTriangular\n";
                    return std::make_shared<CartesianR2ZoniShiftedTriangular>();
                default:
                    throw std::runtime_error("Wrong configuration for the geometry\n");
                }
            case problem_type::POLAR_R6:
                switch (geometry) {
                case geometry_type::CIRCULAR:
                    std::cout << "Using PolarR6ZoniShiftedCircular\n";
                    return std::make_shared<PolarR6ZoniShiftedCircular>();
                case geometry_type::SHAFRANOV:
                    std::cout << "Using PolarR6ZoniShiftedShafranov\n";
                    return std::make_shared<PolarR6ZoniShiftedShafranov>();
                case geometry_type::TRIANGULAR:
                    std::cout << "Using PolarR6ZoniShiftedTriangular\n";
                    return std::make_shared<PolarR6ZoniShiftedTriangular>();
                default:
                    throw std::runtime_error("Wrong configuration for the geometry\n");
                }
            case problem_type::CARTESIAN_R6:
                switch (geometry) {
                case geometry_type::CIRCULAR:
                    std::cout << "Using CartesianR6ZoniShiftedCircular\n";
                    return std::make_shared<CartesianR6ZoniShiftedCircular>();
                case geometry_type::SHAFRANOV:
                    std::cout << "Using CartesianR6ZoniShiftedShafranov\n";
                    return std::make_shared<CartesianR6ZoniShiftedShafranov>();
                case geometry_type::TRIANGULAR:
                    std::cout << "Using CartesianR6ZoniShiftedTriangular\n";
                    return std::make_shared<CartesianR6ZoniShiftedTriangular>();
                default:
                    throw std::runtime_error("Wrong configuration for the geometry\n");
                }
            default:
                throw std::runtime_error("Wrong configuration for the problem\n");
            }
        case alpha_coeff::POISSON:
            switch (problem) {
            case problem_type::CARTESIAN_R2:
                switch (geometry) {
                case geometry_type::CIRCULAR:
                    std::cout << "Using CartesianR2PoissonCircular\n";
                    return std::make_shared<CartesianR2PoissonCircular>();
                case geometry_type::SHAFRANOV:
                    std::cout << "Using CartesianR2PoissonShafranov\n";
                    return std::make_shared<CartesianR2PoissonShafranov>();
                case geometry_type::TRIANGULAR:
                    std::cout << "Using CartesianR2PoissonTriangular\n";
                    return std::make_shared<CartesianR2PoissonTriangular>();
                default:
                    throw std::runtime_error("Wrong configuration for the geometry\n");
                }
            case problem_type::POLAR_R6:
                switch (geometry) {
                case geometry_type::CIRCULAR:
                    std::cout << "Using PolarR6PoissonCircular\n";
                    return std::make_shared<PolarR6PoissonCircular>();
                case geometry_type::SHAFRANOV:
                    std::cout << "Using PolarR6PoissonShafranov\n";
                    return std::make_shared<PolarR6PoissonShafranov>();
                case geometry_type::TRIANGULAR:
                    std::cout << "Using PolarR6PoissonTriangular\n";
                    return std::make_shared<PolarR6PoissonTriangular>();
                default:
                    throw std::runtime_error("Wrong configuration for the geometry\n");
                }
            case problem_type::CARTESIAN_R6:
                switch (geometry) {
                case geometry_type::CIRCULAR:
                    std::cout << "Using CartesianR6PoissonCircular\n";
                    return std::make_shared<CartesianR6PoissonCircular>();
                case geometry_type::SHAFRANOV:
                    std::cout << "Using CartesianR6PoissonShafranov\n";
                    return std::make_shared<CartesianR6PoissonShafranov>();
                case geometry_type::TRIANGULAR:
                    std::cout << "Using CartesianR6PoissonTriangular\n";
                    return std::make_shared<CartesianR6PoissonTriangular>();
                default:
                    throw std::runtime_error("Wrong configuration for the geometry\n");
                }
            default:
                throw std::runtime_error("Wrong configuration for the problem\n");
            }
        default:
            throw std::runtime_error("Wrong configuration for the alpha coefficient\n");
        }
    case beta_coeff::ALPHA_INVERSE:
        switch (alpha) {
        case alpha_coeff::SONNENDRUCKER:
            switch (problem) {
            case problem_type::CARTESIAN_R2:
                switch (geometry) {
                case geometry_type::CIRCULAR:
                    std::cout << "Using CartesianR2GyroSonnendruckerCircular\n";
                    return std::make_shared<CartesianR2GyroSonnendruckerCircular>();
                case geometry_type::SHAFRANOV:
                    std::cout << "Using CartesianR2GyroSonnendruckerShafranov\n";
                    return std::make_shared<CartesianR2GyroSonnendruckerShafranov>();
                case geometry_type::TRIANGULAR:
                    std::cout << "Using CartesianR2GyroSonnendruckerTriangular\n";
                    return std::make_shared<CartesianR2GyroSonnendruckerTriangular>();
                default:
                    throw std::runtime_error("Wrong configuration for the geometry\n");
                }
            case problem_type::POLAR_R6:
                switch (geometry) {
                case geometry_type::CIRCULAR:
                    std::cout << "Using PolarR6GyroSonnendruckerCircular\n";
                    return std::make_shared<PolarR6GyroSonnendruckerCircular>();
                case geometry_type::SHAFRANOV:
                    std::cout << "Using PolarR6GyroSonnendruckerShafranov\n";
                    return std::make_shared<PolarR6GyroSonnendruckerShafranov>();
                case geometry_type::TRIANGULAR:
                    std::cout << "Using PolarR6GyroSonnendruckerTriangular\n";
                    return std::make_shared<PolarR6GyroSonnendruckerTriangular>();
                default:
                    throw std::runtime_error("Wrong configuration for the geometry\n");
                }      
            case problem_type::CARTESIAN_R6:
                switch (geometry) {
                case geometry_type::CIRCULAR:
                    std::cout << "Using CartesianR6GyroSonnendruckerCircular\n";
                    return std::make_shared<CartesianR6GyroSonnendruckerCircular>();
                case geometry_type::SHAFRANOV:
                    std::cout << "Using CartesianR6GyroSonnendruckerShafranov\n";
                    return std::make_shared<CartesianR6GyroSonnendruckerShafranov>();
                case geometry_type::TRIANGULAR:
                    std::cout << "Using CartesianR6GyroSonnendruckerTriangular\n";
                    return std::make_shared<CartesianR6GyroSonnendruckerTriangular>();
                default:
                    throw std::runtime_error("Wrong configuration for the geometry\n"); 
                } 
            default:
                throw std::runtime_error("Wrong configuration for the problem\n");
            }
        case alpha_coeff::ZONI:
            switch (problem) {
            case problem_type::CARTESIAN_R2:
                switch (geometry) {
                case geometry_type::CIRCULAR:
                    std::cout << "Using CartesianR2GyroZoniCircular\n";
                    return std::make_shared<CartesianR2GyroZoniCircular>();
                case geometry_type::SHAFRANOV:
                    std::cout << "Using CartesianR2GyroZoniShafranov\n";
                    return std::make_shared<CartesianR2GyroZoniShafranov>();
                case geometry_type::TRIANGULAR:
                    std::cout << "Using CartesianR2GyroZoniTriangular\n";
                    return std::make_shared<CartesianR2GyroZoniTriangular>();
                default:
                    throw std::runtime_error("Wrong configuration for the geometry\n");
                }
            case problem_type::POLAR_R6:
                switch (geometry) {
                case geometry_type::CIRCULAR:
                    std::cout << "Using PolarR6GyroZoniCircular\n";
                    return std::make_shared<PolarR6GyroZoniCircular>(); 
                case geometry_type::SHAFRANOV:
                    std::cout << "Using PolarR6GyroZoniShafranov\n";
                    return std::make_shared<PolarR6GyroZoniShafranov>();  
                case geometry_type::TRIANGULAR:
                    std::cout << "Using PolarR6GyroZoniTriangular\n";
                    return std::make_shared<PolarR6GyroZoniTriangular>();  
                default:
                    throw std::runtime_error("Wrong configuration for the geometry\n");
                }
            case problem_type::CARTESIAN_R6:
                switch (geometry) {
                case geometry_type::CIRCULAR:
                    std::cout << "Using CartesianR6GyroZoniCircular\n";
                    return std::make_shared<CartesianR6GyroZoniCircular>();
                case geometry_type::SHAFRANOV:
                    std::cout << "Using CartesianR6GyroZoniShafranov\n";
                    return std::make_shared<CartesianR6GyroZoniShafranov>();
                case geometry_type::TRIANGULAR:
                    std::cout << "Using CartesianR6GyroZoniTriangular\n";
                    return std::make_shared<CartesianR6GyroZoniTriangular>();
                default:
                    throw std::runtime_error("Wrong configuration for the geometry\n");
                }
            default:
                throw std::runtime_error("Wrong configuration for the problem\n");  
            }
        case alpha_coeff::ZONI_SHIFTED:
            switch (problem) {
            case problem_type::CARTESIAN_R2:
                switch (geometry) {
                case geometry_type::CIRCULAR:
                    std::cout << "Using CartesianR2GyroZoniShiftedCircular\n";
                    return std::make_shared<CartesianR2GyroZoniShiftedCircular>();
                case geometry_type::SHAFRANOV:
                    std::cout << "Using CartesianR2GyroZoniShiftedShafranov\n";
                    return std::make_shared<CartesianR2GyroZoniShiftedShafranov>();  
                case geometry_type::TRIANGULAR:
                    std::cout << "Using CartesianR2GyroZoniShiftedTriangular\n";
                    return std::make_shared<CartesianR2GyroZoniShiftedTriangular>();      
                default:
                    throw std::runtime_error("Wrong configuration for the geometry\n");         
                }
            case problem_type::POLAR_R6:
                switch (geometry) {
                case geometry_type::CIRCULAR:
                    std::cout << "Using PolarR6GyroZoniShiftedCircular\n";
                    return std::make_shared<PolarR6GyroZoniShiftedCircular>();       
                case geometry_type::SHAFRANOV:
                    std::cout << "Using PolarR6GyroZoniShiftedShafranov\n";
                    return std::make_shared<PolarR6GyroZoniShiftedShafranov>();      
                case geometry_type::TRIANGULAR:
                    std::cout << "Using PolarR6GyroZoniShiftedTriangular\n";
                    return std::make_shared<PolarR6GyroZoniShiftedTriangular>();   
                case geometry_type::CULHAM:
                    std::cout << "Using PolarR6GyroZoniShiftedCulham\n";
                    return std::make_shared<PolarR6GyroZoniShiftedCulham>();
                default:
                    throw std::runtime_error("Wrong configuration for the geometry\n"); 
                }
            case problem_type::CARTESIAN_R6:
                switch (geometry) {
                case geometry_type::CIRCULAR:
                    std::cout << "Using CartesianR6GyroZoniShiftedCircular\n";
                    return std::make_shared<CartesianR6GyroZoniShiftedCircular>();
                case geometry_type::SHAFRANOV:
                    std::cout << "Using CartesianR6GyroZoniShiftedShafranov\n";
                    return std::make_shared<CartesianR6GyroZoniShiftedShafranov>();    
                case geometry_type::TRIANGULAR:
                    std::cout << "Using CartesianR6GyroZoniShiftedTriangular\n";
                    return std::make_shared<CartesianR6GyroZoniShiftedTriangular>();           
                default:
                    throw std::runtime_error("Wrong configuration for the geometry\n");
                } 
            case problem_type::REFINED_RADIUS:
                switch (geometry) {
                case geometry_type::CIRCULAR:
                    std::cout << "Using RefinedGyroZoniShiftedCircular\n";
                    return std::make_shared<RefinedGyroZoniShiftedCircular>();
                case geometry_type::SHAFRANOV:
                    std::cout << "Using RefinedGyroZoniShiftedShafranov\n";
                    return std::make_shared<RefinedGyroZoniShiftedShafranov>();
                case geometry_type::TRIANGULAR:
                    std::cout << "Using RefinedGyroZoniShiftedTriangular\n";
                    return std::make_shared<RefinedGyroZoniShiftedTriangular>();
                case geometry_type::CULHAM:
                    std::cout << "Using RefinedGyroZoniShiftedCulham\n";
                    return std::make_shared<RefinedGyroZoniShiftedCulham>();
                default:
                    throw std::runtime_error("Wrong configuration for the geometry\n");
                }
            default:
                throw std::runtime_error("Wrong configuration for the problem\n");
            }
        case alpha_coeff::POISSON:
            throw std::runtime_error("Beta coeff cannot be 1/0");
        default:
            throw std::runtime_error("Wrong configuration for the alpha coefficient\n");
        }
    default:
        throw std::runtime_error("Wrong configuration for the beta coefficient\n");
    }
    throw std::runtime_error("Error in GMGPolar::selectExactFunctionsClass()\n");
    return nullptr;
}