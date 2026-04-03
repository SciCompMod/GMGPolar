#pragma once

template <class LevelCacheType>
ResidualGive<LevelCacheType>::ResidualGive(const PolarGrid& grid, const LevelCacheType& level_cache,
                                           bool DirBC_Interior, int num_omp_threads)
    : Residual<LevelCacheType>(grid, level_cache, DirBC_Interior, num_omp_threads)
{
}

/* ------------ */
/* result = A*x */
/* ------------ */
template <class LevelCacheType>
void ResidualGive<LevelCacheType>::applySystemOperator(Vector<double> result, ConstVector<double> x) const
{
    assert(result.size() == x.size());

    assign(result, 0.0);

    const PolarGrid& grid     = Residual<LevelCacheType>::grid_;
    const int num_omp_threads = Residual<LevelCacheType>::num_omp_threads_;

    const int num_smoother_circles    = grid.numberSmootherCircles();
    const int additional_radial_tasks = grid.ntheta() % 3;
    const int num_radial_tasks        = grid.ntheta() - additional_radial_tasks;

    /* ---------------- */
    /* Circular section */
    /* ---------------- */
    // We parallelize the loop with step 3 to avoid data race conditions between adjacent circles.
#pragma omp parallel num_threads(num_omp_threads)
    {
#pragma omp for
        for (int i_r = 0; i_r < num_smoother_circles; i_r += 3) {
            applyCircleSection(i_r, result, x);
        } /* Implicit barrier */
#pragma omp for
        for (int i_r = 1; i_r < num_smoother_circles; i_r += 3) {
            applyCircleSection(i_r, result, x);
        } /* Implicit barrier */
#pragma omp for
        for (int i_r = 2; i_r < num_smoother_circles; i_r += 3) {
            applyCircleSection(i_r, result, x);
        } /* Implicit barrier */
    }

    /* ---------------- */
    /* Radial section */
    /* ---------------- */
    // We parallelize the loop with step 3 to avoid data race conditions between adjacent radial lines.
    // Due to the periodicity in the angular direction, we can have at most 2 additional radial tasks
    // that are handled serially before the parallel loops.
    if (additional_radial_tasks > 0) {
        const int i_theta = 0;
        applyRadialSection(i_theta, result, x);
    }

    if (additional_radial_tasks > 1) {
        const int i_theta = 1;
        applyRadialSection(i_theta, result, x);
    }

#pragma omp parallel num_threads(num_omp_threads)
    {
#pragma omp for
        for (int radial_task = 0; radial_task < num_radial_tasks; radial_task += 3) {
            const int i_theta = radial_task + additional_radial_tasks;
            applyRadialSection(i_theta, result, x);
        } /* Implicit barrier */
#pragma omp for
        for (int radial_task = 1; radial_task < num_radial_tasks; radial_task += 3) {
            const int i_theta = radial_task + additional_radial_tasks;
            applyRadialSection(i_theta, result, x);
        } /* Implicit barrier */
#pragma omp for
        for (int radial_task = 2; radial_task < num_radial_tasks; radial_task += 3) {
            const int i_theta = radial_task + additional_radial_tasks;
            applyRadialSection(i_theta, result, x);
        } /* Implicit barrier */
    }
}

/* ------------------ */
/* result = rhs - A*x */
/* ------------------ */
template <class LevelCacheType>
void ResidualGive<LevelCacheType>::computeResidual(Vector<double> result, ConstVector<double> rhs,
                                                   ConstVector<double> x) const
{
    assert(result.size() == x.size());

    const int num_omp_threads = Residual<LevelCacheType>::num_omp_threads_;

    applySystemOperator(result, x);

    // Subtract A*x from rhs to get the residual.
    const int n = result.size();
#pragma omp parallel for num_threads(num_omp_threads)
    for (int i = 0; i < n; i++) {
        result[i] = rhs[i] - result[i];
    }
}
