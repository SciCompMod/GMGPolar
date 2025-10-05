#include "../../../include/SystemOperator/SystemOperatorGive/system_operator.h"
#include "../../../include/SystemOperator/system_operator.h"

SystemOperatorGive::SystemOperatorGive(const PolarGrid& grid, const LevelCache& level_cache,
                                       const DomainGeometry& domain_geometry,
                                       const DensityProfileCoefficients& density_profile_coefficients,
                                       bool DirBC_Interior, int num_omp_threads)
    : SystemOperator(grid, level_cache, domain_geometry, density_profile_coefficients, DirBC_Interior, num_omp_threads)
{
}

/* ------------------ */
/* result = rhs - A*x */

// clang-format off
void SystemOperatorGive::apply(Vector<double>& result, const Vector<double>& x) const
{
    assert(result.size() == x.size());

    assign(result, 0.0);

    if (num_omp_threads_ == 1) {
        /* Single-threaded execution */
        for (int i_r = 0; i_r < grid_.numberSmootherCircles(); i_r++) {
            applyCircleSection(i_r, result, x);
        }
        for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++) {
            applyRadialSection(i_theta, result, x);
        }
    }
    else {
        /* Multi-threaded execution */
        const int num_circle_tasks        = grid_.numberSmootherCircles();
        const int additional_radial_tasks = grid_.ntheta() % 3;
        const int num_radial_tasks        = grid_.ntheta() - additional_radial_tasks;

        #pragma omp parallel num_threads(num_omp_threads_)
        {
            /* Circle Section 0 */
            #pragma omp for
            for (int circle_task = 0; circle_task < num_circle_tasks; circle_task += 3) {
                int i_r = grid_.numberSmootherCircles() - circle_task - 1;
                applyCircleSection(i_r, result, x);
            }
            /* Circle Section 1 */
            #pragma omp for
            for (int circle_task = 1; circle_task < num_circle_tasks; circle_task += 3) {
                int i_r = grid_.numberSmootherCircles() - circle_task - 1;
                applyCircleSection(i_r, result, x);
            }
            /* Circle Section 2 */
            #pragma omp for nowait
            for (int circle_task = 2; circle_task < num_circle_tasks; circle_task += 3) {
                int i_r = grid_.numberSmootherCircles() - circle_task - 1;
                applyCircleSection(i_r, result, x);
            }

            /* Radial Section 0 */
            #pragma omp for
            for (int radial_task = 0; radial_task < num_radial_tasks; radial_task += 3) {
                if (radial_task > 0) {
                    int i_theta = radial_task + additional_radial_tasks;
                    applyRadialSection(i_theta, result, x);
                }
                else {
                    if (additional_radial_tasks == 0) {
                        applyRadialSection(0, result, x);
                    }
                    else if (additional_radial_tasks >= 1) {
                        applyRadialSection(0, result, x);
                        applyRadialSection(1, result, x);
                    }
                }
            }
            /* Radial Section 1 */
            #pragma omp for
            for (int radial_task = 1; radial_task < num_radial_tasks; radial_task += 3) {
                if (radial_task > 1) {
                    int i_theta = radial_task + additional_radial_tasks;
                    applyRadialSection(i_theta, result, x);
                }
                else {
                    if (additional_radial_tasks == 0) {
                        applyRadialSection(1, result, x);
                    }
                    else if (additional_radial_tasks == 1) {
                        applyRadialSection(2, result, x);
                    }
                    else if (additional_radial_tasks == 2) {
                        applyRadialSection(2, result, x);
                        applyRadialSection(3, result, x);
                    }
                }
            }
            /* Radial Section 2 */
            #pragma omp for
            for (int radial_task = 2; radial_task < num_radial_tasks; radial_task += 3) {
                int i_theta = radial_task + additional_radial_tasks;
                applyRadialSection(i_theta, result, x);
            }
        }
    }
}
// clang-format on
