#include "../../include/Interpolation/interpolation.h"

#define FINE_NODE_FMG_INTERPOLATION() \
do { \
    if(i_r == 0 || i_r == fineGrid.nr() - 1){ \
        if(i_theta & 1){ \
            double k0 = coarseGrid.angularSpacing(i_theta_coarse-1); \
            double k1 = fineGrid.angularSpacing(i_theta-1); \
            double k2 = fineGrid.angularSpacing(i_theta); \
            double k3 = coarseGrid.angularSpacing(i_theta_coarse+1); \
            \
                double w_theta0 = - k1/k0 * k2/(k0+k1+k2) * (k2+k3)/(k0+k1+k2+k3); \
                double w_theta1 = (k0+k1)/k0 * k2/(k1+k2) * (k2+k3)/(k1+k2+k3); \
                double w_theta2 = (k0+k1)/(k0+k1+k2) * k1/(k1+k2) * (k2+k3)/k3; \
                double w_theta3 = - (k0+k1)/(k0+k1+k2+k3) * k1/(k1+k2+k3) * k2/k3; \
            \
            result[fineGrid.index(i_r, i_theta)] = ( \
                w_theta0 * x[coarseGrid.index(i_r_coarse, i_theta_coarse-1)] + /* (0, -3) */ \
                w_theta1 * x[coarseGrid.index(i_r_coarse, i_theta_coarse  )] + /* (0, -1) */ \
                w_theta2 * x[coarseGrid.index(i_r_coarse, i_theta_coarse+1)] + /* (0, +1) */ \
                w_theta3 * x[coarseGrid.index(i_r_coarse, i_theta_coarse+2)] /* (0, +3) */ \
            ); \
        } \
        else{ \
            result[fineGrid.index(i_r, i_theta)] = \
                x[coarseGrid.index(i_r_coarse, i_theta_coarse)]; /* center */ \
        } \
    } \
    else if(i_r == 1 || i_r == fineGrid.nr() - 2){ \
        if(i_theta & 1){ \
            double k0 = coarseGrid.angularSpacing(i_theta_coarse-1); \
            double k1 = fineGrid.angularSpacing(i_theta-1); \
            double k2 = fineGrid.angularSpacing(i_theta); \
            double k3 = coarseGrid.angularSpacing(i_theta_coarse+1); \
            \
                double w_theta0 = - k1/k0 * k2/(k0+k1+k2) * (k2+k3)/(k0+k1+k2+k3); \
                double w_theta1 = (k0+k1)/k0 * k2/(k1+k2) * (k2+k3)/(k1+k2+k3); \
                double w_theta2 = (k0+k1)/(k0+k1+k2) * k1/(k1+k2) * (k2+k3)/k3; \
                double w_theta3 = - (k0+k1)/(k0+k1+k2+k3) * k1/(k1+k2+k3) * k2/k3; \
            \
            double left_value = ( \
                w_theta0 * x[coarseGrid.index(i_r_coarse, i_theta_coarse-1)] + /* (-1, -3) */ \
                w_theta1 * x[coarseGrid.index(i_r_coarse, i_theta_coarse  )] + /* (-1, -1) */ \
                w_theta2 * x[coarseGrid.index(i_r_coarse, i_theta_coarse+1)] + /* (-1, +1) */ \
                w_theta3 * x[coarseGrid.index(i_r_coarse, i_theta_coarse+2)]   /* (-1, +3) */ \
            ); \
            double right_value = ( \
                w_theta0 * x[coarseGrid.index(i_r_coarse+1, i_theta_coarse-1)] + /* (+1, -3) */ \
                w_theta1 * x[coarseGrid.index(i_r_coarse+1, i_theta_coarse  )] + /* (+1, -1) */ \
                w_theta2 * x[coarseGrid.index(i_r_coarse+1, i_theta_coarse+1)] + /* (+1, +1) */ \
                w_theta3 * x[coarseGrid.index(i_r_coarse+1, i_theta_coarse+2)]   /* (+1, +3) */ \
            ); \
            \
            double h1 = fineGrid.radialSpacing(i_r-1); \
            double h2 = fineGrid.radialSpacing(i_r); \
            result[fineGrid.index(i_r, i_theta)] = (h1 * left_value + h2 * right_value) / (h1 + h2); \
        } \
        else{ \
            double h1 = fineGrid.radialSpacing(i_r-1); \
            double h2 = fineGrid.radialSpacing(i_r); \
            result[fineGrid.index(i_r, i_theta)] = ( \
                h1 * x[coarseGrid.index(i_r_coarse, i_theta_coarse)] + /* left */ \
                h2 * x[coarseGrid.index(i_r_coarse+1, i_theta_coarse)] /* right */ \
            ) / (h1 + h2); \
        } \
    } \
    else{ \
        if(i_r & 1){ \
            if(i_theta & 1){ \
                double k0 = coarseGrid.angularSpacing(i_theta_coarse-1); \
                double k1 = fineGrid.angularSpacing(i_theta-1); \
                double k2 = fineGrid.angularSpacing(i_theta); \
                double k3 = coarseGrid.angularSpacing(i_theta_coarse+1); \
                \
                double w_theta0 = - k1/k0 * k2/(k0+k1+k2) * (k2+k3)/(k0+k1+k2+k3); \
                double w_theta1 = (k0+k1)/k0 * k2/(k1+k2) * (k2+k3)/(k1+k2+k3); \
                double w_theta2 = (k0+k1)/(k0+k1+k2) * k1/(k1+k2) * (k2+k3)/k3; \
                double w_theta3 = - (k0+k1)/(k0+k1+k2+k3) * k1/(k1+k2+k3) * k2/k3; \
                \
                double outer_left_value = ( \
                    w_theta0 * x[coarseGrid.index(i_r_coarse-1, i_theta_coarse-1)] + /* (-3, -3) */ \
                    w_theta1 * x[coarseGrid.index(i_r_coarse-1, i_theta_coarse  )] + /* (-3, -1) */ \
                    w_theta2 * x[coarseGrid.index(i_r_coarse-1, i_theta_coarse+1)] + /* (-3, +1) */ \
                    w_theta3 * x[coarseGrid.index(i_r_coarse-1, i_theta_coarse+2)]   /* (-3, +3) */ \
                ); \
                double inner_left_value = ( \
                    w_theta0 * x[coarseGrid.index(i_r_coarse, i_theta_coarse-1)] + /* (-1, -3) */ \
                    w_theta1 * x[coarseGrid.index(i_r_coarse, i_theta_coarse  )] + /* (-1, -1) */ \
                    w_theta2 * x[coarseGrid.index(i_r_coarse, i_theta_coarse+1)] + /* (-1, +1) */ \
                    w_theta3 * x[coarseGrid.index(i_r_coarse, i_theta_coarse+2)]   /* (-1, +3) */ \
                ); \
                double inner_right_value = ( \
                    w_theta0 * x[coarseGrid.index(i_r_coarse+1, i_theta_coarse-1)] + /* (+1, -3) */ \
                    w_theta1 * x[coarseGrid.index(i_r_coarse+1, i_theta_coarse  )] + /* (+1, -1) */ \
                    w_theta2 * x[coarseGrid.index(i_r_coarse+1, i_theta_coarse+1)] + /* (+1, +1) */ \
                    w_theta3 * x[coarseGrid.index(i_r_coarse+1, i_theta_coarse+2)]   /* (+1, +3) */ \
                ); \
                double outer_right_value = ( \
                    w_theta0 * x[coarseGrid.index(i_r_coarse+2, i_theta_coarse-1)] + /* (+3, -3) */ \
                    w_theta1 * x[coarseGrid.index(i_r_coarse+2, i_theta_coarse  )] + /* (+3, -1) */ \
                    w_theta2 * x[coarseGrid.index(i_r_coarse+2, i_theta_coarse+1)] + /* (+3, +1) */ \
                    w_theta3 * x[coarseGrid.index(i_r_coarse+2, i_theta_coarse+2)]   /* (+3, +3) */ \
                ); \
                \
                double h0 = coarseGrid.radialSpacing(i_r_coarse-1); \
                double h1 = fineGrid.radialSpacing(i_r-1); \
                double h2 = fineGrid.radialSpacing(i_r); \
                double h3 = coarseGrid.radialSpacing(i_r_coarse+1); \
                \
                double w_r0 = - h1/h0 * h2/(h0+h1+h2) * (h2+h3)/(h0+h1+h2+h3); \
                double w_r1 = (h0+h1)/h0 * h2/(h1+h2) * (h2+h3)/(h1+h2+h3); \
                double w_r2 = (h0+h1)/(h0+h1+h2) * h1/(h1+h2) * (h2+h3)/h3; \
                double w_r3 = - (h0+h1)/(h0+h1+h2+h3) * h1/(h1+h2+h3) * h2/h3; \
                \
                result[fineGrid.index(i_r, i_theta)] = ( \
                    w_r0 * outer_left_value + \
                    w_r1 * inner_left_value + \
                    w_r2 * inner_right_value + \
                    w_r3 * outer_right_value \
                ); \
            } \
            else{ \
                double h0 = coarseGrid.radialSpacing(i_r_coarse-1); \
                double h1 = fineGrid.radialSpacing(i_r-1); \
                double h2 = fineGrid.radialSpacing(i_r); \
                double h3 = coarseGrid.radialSpacing(i_r_coarse+1); \
                \
                double w_r0 = - h1/h0 * h2/(h0+h1+h2) * (h2+h3)/(h0+h1+h2+h3); \
                double w_r1 = (h0+h1)/h0 * h2/(h1+h2) * (h2+h3)/(h1+h2+h3); \
                double w_r2 = (h0+h1)/(h0+h1+h2) * h1/(h1+h2) * (h2+h3)/h3; \
                double w_r3 = - (h0+h1)/(h0+h1+h2+h3) * h1/(h1+h2+h3) * h2/h3; \
                \
                result[fineGrid.index(i_r, i_theta)] = ( \
                    w_r0 * x[coarseGrid.index(i_r_coarse-1, i_theta_coarse)] + /* (-3, 0) */ \
                    w_r1 * x[coarseGrid.index(i_r_coarse, i_theta_coarse)] + /* (-1, 0) */ \
                    w_r2 * x[coarseGrid.index(i_r_coarse+1, i_theta_coarse)] + /* (+1, 0) */ \
                    w_r3 * x[coarseGrid.index(i_r_coarse+2, i_theta_coarse)] /* (+3, 0) */ \
                ); \
            } \
        } \
        else{ \
            if(i_theta & 1){ \
                double k0 = coarseGrid.angularSpacing(i_theta_coarse-1); \
                double k1 = fineGrid.angularSpacing(i_theta-1); \
                double k2 = fineGrid.angularSpacing(i_theta); \
                double k3 = coarseGrid.angularSpacing(i_theta_coarse+1); \
                \
                double w_theta0 = - k1/k0 * k2/(k0+k1+k2) * (k2+k3)/(k0+k1+k2+k3); \
                double w_theta1 = (k0+k1)/k0 * k2/(k1+k2) * (k2+k3)/(k1+k2+k3); \
                double w_theta2 = (k0+k1)/(k0+k1+k2) * k1/(k1+k2) * (k2+k3)/k3; \
                double w_theta3 = - (k0+k1)/(k0+k1+k2+k3) * k1/(k1+k2+k3) * k2/k3; \
                \
                result[fineGrid.index(i_r, i_theta)]  = ( \
                    w_theta0 * x[coarseGrid.index(i_r_coarse, i_theta_coarse-1)] + /* (0, -3) */ \
                    w_theta1 * x[coarseGrid.index(i_r_coarse, i_theta_coarse  )] + /* (0, -1) */ \
                    w_theta2 * x[coarseGrid.index(i_r_coarse, i_theta_coarse+1)] + /* (0, +1) */ \
                    w_theta3 * x[coarseGrid.index(i_r_coarse, i_theta_coarse+2)]   /* (0, +3) */ \
                ); \
            } \
            else{ \
                result[fineGrid.index(i_r, i_theta)] = x[coarseGrid.index(i_r_coarse, i_theta_coarse)]; /* center */ \
            } \
        } \
    } \
} while(0)



void Interpolation::applyFMGInterpolation(const Level& fromLevel, const Level& toLevel, Vector<double>& result, const Vector<double>& x) const{
    assert(toLevel.level() == fromLevel.level() - 1);

    omp_set_num_threads(threads_per_level_[toLevel.level()]);

    const PolarGrid& coarseGrid = fromLevel.grid();
    const PolarGrid& fineGrid = toLevel.grid();

    assert(x.size() == coarseGrid.numberOfNodes());
    assert(result.size() == fineGrid.numberOfNodes());

    #pragma omp parallel if(fineGrid.numberOfNodes() > 10'000)
    {
        /* Circluar Indexing Section */
        /* For loop matches circular access pattern */
        #pragma omp for nowait
        for (int i_r = 0; i_r < fineGrid.numberSmootherCircles(); i_r++){
            int i_r_coarse = i_r >> 1;
            for (int i_theta = 0; i_theta < fineGrid.ntheta(); i_theta++){
                int i_theta_coarse = i_theta >> 1;
                FINE_NODE_FMG_INTERPOLATION();
            }
        }

        /* Radial Indexing Section */
        /* For loop matches radial access pattern */
        #pragma omp for nowait
        for (int i_theta = 0; i_theta < fineGrid.ntheta(); i_theta++){
            int i_theta_coarse = i_theta >> 1;
            for (int i_r = fineGrid.numberSmootherCircles(); i_r < fineGrid.nr(); i_r++){
                int i_r_coarse = i_r >> 1;
                FINE_NODE_FMG_INTERPOLATION();
            }
        }
    }
}