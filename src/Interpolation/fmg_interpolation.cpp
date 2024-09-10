#include "../../include/Interpolation/interpolation.h"

#define FINE_NODE_PROLONGATION() \
do { \
    if(i_r == 0 || i_r == fineGrid.nr() - 1){ \
        if(i_theta & 1){ \
            result[fineGrid.index(i_r, i_theta)] = 0.5 * ( \
                x[coarseGrid.index(i_r_coarse, i_theta_coarse)] + /* bottom */ \
                x[coarseGrid.index(i_r_coarse, i_theta_coarse+1)] /* top */ \
            ); \
        } \
        else{ \
            result[fineGrid.index(i_r, i_theta)] = \
                x[coarseGrid.index(i_r_coarse, i_theta_coarse)]; /* center */ \
        } \
    } \
    else if(i_r == 1 || i_r == fineGrid.nr() - 2){ \
        if(i_theta & 1){ \
            double left_value = 1.0 / 16.0 * ( \
                -1.0 * x[coarseGrid.index(i_r_coarse, i_theta_coarse-1)] + /* (-1, -3) */ \
                9.0 * x[coarseGrid.index(i_r_coarse, i_theta_coarse)] + /* (-1, -1) */ \
                9.0 * x[coarseGrid.index(i_r_coarse, i_theta_coarse+1)] + /* (-1, +1) */ \
                -1.0 * x[coarseGrid.index(i_r_coarse, i_theta_coarse+2)] /* (-1, +3) */ \
            ); \
            double right_value = 1.0 / 16.0 * ( \
                -1.0 * + x[coarseGrid.index(i_r_coarse+1, i_theta_coarse-1)] + /* (+1, -3) */ \
                9.0 * + x[coarseGrid.index(i_r_coarse+1, i_theta_coarse)] + /* (+1, -1) */ \
                9.0 * + x[coarseGrid.index(i_r_coarse+1, i_theta_coarse+1)] + /* (+1, +1) */ \
                -1.0 * x[coarseGrid.index(i_r_coarse+1, i_theta_coarse+2)] /* (+1, +3) */ \
            ); \
            result[fineGrid.index(i_r, i_theta)] = 0.5 * ( \
                left_value + /* left */ \
                right_value /* right */ \
            ); \
        } \
        else{ \
            result[fineGrid.index(i_r, i_theta)] = 0.5 * ( \
                x[coarseGrid.index(i_r_coarse, i_theta_coarse)] + /* left */ \
                x[coarseGrid.index(i_r_coarse+1, i_theta_coarse)] /* right */ \
            ); \
        } \
    } \
    else{ \
        if(i_r & 1){ \
            if(i_theta & 1){ \
                double left_first_value = 1.0 / 16.0 * ( \
                    -1.0 * x[coarseGrid.index(i_r_coarse-1, i_theta_coarse-1)] + /* (-3, -3) */ \
                    9.0 * x[coarseGrid.index(i_r_coarse-1, i_theta_coarse)] + /* (-3, -1) */ \
                    9.0 * x[coarseGrid.index(i_r_coarse-1, i_theta_coarse+1)] + /* (-3, +1) */ \
                    -1.0 * x[coarseGrid.index(i_r_coarse-1, i_theta_coarse+2)] /* (-3, +3) */ \
                ); \
                double left_second_value = 1.0 / 16.0 * ( \
                    -1.0 * x[coarseGrid.index(i_r_coarse, i_theta_coarse-1)] + /* (-1, -3) */ \
                    9.0 * x[coarseGrid.index(i_r_coarse, i_theta_coarse)] + /* (-1, -1) */ \
                    9.0 * x[coarseGrid.index(i_r_coarse, i_theta_coarse+1)] + /* (-1, +1) */ \
                    -1.0 * x[coarseGrid.index(i_r_coarse, i_theta_coarse+2)] /* (-1, +3) */ \
                ); \
                double right_first_value = 1.0 / 16.0 * ( \
                    -1.0 * + x[coarseGrid.index(i_r_coarse+1, i_theta_coarse-1)] + /* (+1, -3) */ \
                    9.0 * + x[coarseGrid.index(i_r_coarse+1, i_theta_coarse)] + /* (+1, -1) */ \
                    9.0 * + x[coarseGrid.index(i_r_coarse+1, i_theta_coarse+1)] + /* (+1, +1) */ \
                    -1.0 * x[coarseGrid.index(i_r_coarse+1, i_theta_coarse+2)] /* (+1, +3) */ \
                ); \
                double right_second_value = 1.0 / 16.0 * ( \
                    -1.0 * + x[coarseGrid.index(i_r_coarse+2, i_theta_coarse-1)] + /* (+3, -3) */ \
                    9.0 * + x[coarseGrid.index(i_r_coarse+2, i_theta_coarse)] + /* (+3, -1) */ \
                    9.0 * + x[coarseGrid.index(i_r_coarse+2, i_theta_coarse+1)] + /* (+3, +1) */ \
                    -1.0 * x[coarseGrid.index(i_r_coarse+2, i_theta_coarse+2)] /* (+3, +3) */ \
                ); \
                result[fineGrid.index(i_r, i_theta)] = 1.0 / 16.0 * ( \
                    -1.0 * left_first_value + \
                    9.0 * left_second_value + \
                    9.0 * right_first_value + \
                    -1.0 * right_second_value \
                ); \
            } \
            else{ \
                result[fineGrid.index(i_r, i_theta)] = 1.0 / 16.0 * ( \
                    -1.0 * + x[coarseGrid.index(i_r_coarse-1, i_theta_coarse)] + /* (-3, 0) */ \
                    9.0 * + x[coarseGrid.index(i_r_coarse, i_theta_coarse)] + /* (-1, 0) */ \
                    9.0 * + x[coarseGrid.index(i_r_coarse+1, i_theta_coarse)] + /* (+1, 0) */ \
                    -1.0 * x[coarseGrid.index(i_r_coarse+2, i_theta_coarse)] /* (+3, 0) */ \
                ); \
            } \
        } \
        else{ \
            if(i_theta & 1){ \
                result[fineGrid.index(i_r, i_theta)] = 1.0 / 16.0 * ( \
                    -1.0 * + x[coarseGrid.index(i_r_coarse, i_theta_coarse-1)] + /* (0, -3) */ \
                    9.0 * + x[coarseGrid.index(i_r_coarse, i_theta_coarse)] + /* (0, -1) */ \
                    9.0 * + x[coarseGrid.index(i_r_coarse, i_theta_coarse+1)] + /* (0, +1) */ \
                    -1.0 * x[coarseGrid.index(i_r_coarse, i_theta_coarse+2)] /* (0, +3) */ \
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
                FINE_NODE_PROLONGATION();
            }
        }

        /* Radial Indexing Section */
        /* For loop matches radial access pattern */
        #pragma omp for nowait
        for (int i_theta = 0; i_theta < fineGrid.ntheta(); i_theta++){
            int i_theta_coarse = i_theta >> 1;
            for (int i_r = fineGrid.numberSmootherCircles(); i_r < fineGrid.nr(); i_r++){
                int i_r_coarse = i_r >> 1;
                FINE_NODE_PROLONGATION();
            }
        }
    }
}