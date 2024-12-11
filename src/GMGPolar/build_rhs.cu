#include "../../include/GMGPolar/gmgpolar.h"


void GMGPolar::build_rhs_f(const Level& level, Vector<double>& rhs_f)
{
    const PolarGrid& grid = level.grid();
    assert(rhs_f.size() == grid.numberOfNodes());

    const auto& sin_theta_cache = level.levelCache().sin_theta();
    const auto& cos_theta_cache = level.levelCache().cos_theta();

    #pragma omp parallel
    {
        // ----------------------------------------- //
        // Store rhs values (circular index section) //
        // ----------------------------------------- //
        #pragma omp for nowait
        for (int i_r = 0; i_r < grid.numberSmootherCircles(); i_r++)
        {
            double r = grid.radius(i_r);
            for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++)
            {
                double theta = grid.theta(i_theta);
                double sin_theta = sin_theta_cache[i_theta];
                double cos_theta = cos_theta_cache[i_theta];

                if ((0 < i_r && i_r < grid.nr() - 1) || (i_r == 0 && !DirBC_Interior_))
                {
                    rhs_f[grid.index(i_r, i_theta)] = source_term_->rhs_f(r, theta, sin_theta, cos_theta);
                }
                else if (i_r == 0 && DirBC_Interior_)
                {
                    rhs_f[grid.index(i_r, i_theta)] = boundary_conditions_->u_D_Interior(r, theta, sin_theta, cos_theta);
                }
                else if (i_r == grid.nr() - 1)
                {
                    rhs_f[grid.index(i_r, i_theta)] = boundary_conditions_->u_D(r, theta, sin_theta, cos_theta);
                }
            }
        }

        // --------------------------------------- //
        // Store rhs values (radial index section) //
        // --------------------------------------- //
        #pragma omp for
        for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++)
        {
            double theta = grid.theta(i_theta);
            double sin_theta = sin_theta_cache[i_theta];
            double cos_theta = cos_theta_cache[i_theta];

            for (int i_r = grid.numberSmootherCircles(); i_r < grid.nr(); i_r++)
            {
                double r = grid.radius(i_r);
                if ((0 < i_r && i_r < grid.nr() - 1) || (i_r == 0 && !DirBC_Interior_))
                {
                    rhs_f[grid.index(i_r, i_theta)] = source_term_->rhs_f(r, theta, sin_theta, cos_theta);
                }
                else if (i_r == 0 && DirBC_Interior_)
                {
                    rhs_f[grid.index(i_r, i_theta)] = boundary_conditions_->u_D_Interior(r, theta, sin_theta, cos_theta);
                }
                else if (i_r == grid.nr() - 1)
                {
                    rhs_f[grid.index(i_r, i_theta)] = boundary_conditions_->u_D(r, theta, sin_theta, cos_theta);
                }
            }
        }
    }
}

void GMGPolar::discretize_rhs_f(const Level& level, Vector<double>& rhs_f)
{
    const PolarGrid& grid = level.grid();
    assert(rhs_f.size() == grid.numberOfNodes());

    const auto& detDF_cache = level.levelCache().detDF();
    #pragma omp parallel
    {
        // ---------------------------------------------- //
        // Discretize rhs values (circular index section) //
        // ---------------------------------------------- //
        #pragma omp for nowait
        for (int i_r = 0; i_r < grid.numberSmootherCircles(); i_r++)
        {
            double r = grid.radius(i_r);
            for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++)
            {
                double theta = grid.theta(i_theta);
                if ((0 < i_r && i_r < grid.nr() - 1) || (i_r == 0 && !DirBC_Interior_))
                {
                    double h1 = (i_r == 0) ? 2.0 * grid.radius(0) : grid.radialSpacing(i_r - 1);
                    double h2 = grid.radialSpacing(i_r);
                    double k1 = grid.angularSpacing(i_theta - 1);
                    double k2 = grid.angularSpacing(i_theta);
                    const double detDF = detDF_cache[grid.index(i_r, i_theta)];
                    rhs_f[grid.index(i_r, i_theta)] *= 0.25 * (h1 + h2) * (k1 + k2) * fabs(detDF);
                }
                else if (i_r == 0 && DirBC_Interior_)
                {
                    rhs_f[grid.index(i_r, i_theta)] *= 1.0;
                }
                else if (i_r == grid.nr() - 1)
                {
                    rhs_f[grid.index(i_r, i_theta)] *= 1.0;
                }
            }
        }

        // -------------------------------------------- //
        // Discretize rhs values (radial index section) //
        // -------------------------------------------- //
        #pragma omp for nowait
        for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++)
        {
            double theta = grid.theta(i_theta);
            for (int i_r = grid.numberSmootherCircles(); i_r < grid.nr(); i_r++)
            {
                double r = grid.radius(i_r);
                if ((0 < i_r && i_r < grid.nr() - 1) || (i_r == 0 && !DirBC_Interior_))
                {
                    double h1 = (i_r == 0) ? 2.0 * grid.radius(0) : grid.radialSpacing(i_r - 1);
                    double h2 = grid.radialSpacing(i_r);
                    double k1 = grid.angularSpacing(i_theta - 1);
                    double k2 = grid.angularSpacing(i_theta);
                    const double detDF = detDF_cache[grid.index(i_r, i_theta)];
                    rhs_f[grid.index(i_r, i_theta)] *= 0.25 * (h1 + h2) * (k1 + k2) * fabs(detDF);
                }
                else if (i_r == 0 && DirBC_Interior_)
                {
                    rhs_f[grid.index(i_r, i_theta)] *= 1.0;
                }
                else if (i_r == grid.nr() - 1)
                {
                    rhs_f[grid.index(i_r, i_theta)] *= 1.0;
                }
            }
        }
    }
}



__global__ void build_rhs_f_kernel(PolarGrid* grid, double* rhs_f, double* sin_theta_cache, double* cos_theta_cache, bool DirBC_Interior, SourceTerm* source_term, BoundaryConditions* boundary_conditions) {
    int i_r = blockIdx.x * blockDim.x + threadIdx.x;
    int i_theta = blockIdx.y * blockDim.y + threadIdx.y;

    if (i_r >= grid->nr() || i_theta >= grid->ntheta()) return;

    double r = grid->radius(i_r);
    double theta = grid->theta(i_theta);

    double sin_theta = sin_theta_cache[i_theta];
    double cos_theta = cos_theta_cache[i_theta];

    double value;

    if ((0 < i_r && i_r < grid->nr() - 1) || (i_r == 0 && !DirBC_Interior))
    {
        value = source_term->rhs_f(r, theta, sin_theta, cos_theta);
    }
    else if (i_r == 0 && DirBC_Interior)
    {
        value = boundary_conditions->u_D_Interior(r, theta, sin_theta, cos_theta);
    }
    else if (i_r == grid->nr() - 1)
    {
        value = boundary_conditions->u_D(r, theta, sin_theta, cos_theta);
    }

    rhs_f[grid->index(i_r, i_theta)] = value;
}

void GMGPolar::build_rhs_f(const Level& level, GPU_Vector<double>& rhs_f)
{
    const PolarGrid& grid = level.grid();
    assert(rhs_f.size() == grid.numberOfNodes());

    const GPU_Vector<double>& sin_theta_cache = level.levelCache().GPU_sin_theta();
    const GPU_Vector<double>& cos_theta_cache = level.levelCache().GPU_cos_theta();

    dim3 threadsPerBlock(16, 16);
    dim3 numBlocks((grid.nr() + threadsPerBlock.x - 1) / threadsPerBlock.x,
                   (grid.ntheta() + threadsPerBlock.y - 1) / threadsPerBlock.y);

    SourceTerm* device_source_term;
    cudaMalloc(&device_source_term, sizeof(SourceTerm));
    cudaMemcpy(device_source_term, source_term_.get(), sizeof(SourceTerm), cudaMemcpyHostToDevice);
    BoundaryConditions* device_boundary_conditions;
    cudaMalloc(&device_boundary_conditions, sizeof(BoundaryConditions));
    cudaMemcpy(device_boundary_conditions, boundary_conditions_.get(), sizeof(BoundaryConditions), cudaMemcpyHostToDevice);

    build_rhs_f_kernel<<<numBlocks, threadsPerBlock>>>(
        level.device_grid(), rhs_f.data(), sin_theta_cache.data(), cos_theta_cache.data(), DirBC_Interior_, device_source_term, device_boundary_conditions);

    cudaDeviceSynchronize();

    cudaFree(device_source_term);
    cudaFree(device_boundary_conditions);
}

__global__ void discretize_rhs_f_kernel(PolarGrid* grid, double* rhs_f, double* sin_theta_cache, double* cos_theta_cache, bool DirBC_Interior, DomainGeometry* domain_geometry) {
    int i_r = blockIdx.x * blockDim.x + threadIdx.x;
    int i_theta = blockIdx.y * blockDim.y + threadIdx.y;

    if (i_r >= grid->nr() || i_theta >= grid->ntheta()) return;

    double r = grid->radius(i_r);
    double theta = grid->theta(i_theta);

    double sin_theta = sin_theta_cache[i_theta];
    double cos_theta = cos_theta_cache[i_theta];

    double value;

    if ((0 < i_r && i_r < grid->nr() - 1) || (i_r == 0 && !DirBC_Interior))
    {
        double h1 = (i_r == 0) ? 2.0 * grid->radius(0) : grid->radialSpacing(i_r - 1);
        double h2 = grid->radialSpacing(i_r);
        double k1 = grid->angularSpacing(i_theta - 1);
        double k2 = grid->angularSpacing(i_theta);
        /* Calculate the elements of the Jacobian matrix for the transformation mapping */
        /* The Jacobian matrix is: */
        /* [Jrr, Jrt] */
        /* [Jtr, Jtt] */
        double Jrr = domain_geometry->dFx_dr(r, theta, sin_theta, cos_theta);
        double Jtr = domain_geometry->dFy_dr(r, theta, sin_theta, cos_theta);
        double Jrt = domain_geometry->dFx_dt(r, theta, sin_theta, cos_theta);
        double Jtt = domain_geometry->dFy_dt(r, theta, sin_theta, cos_theta);
        /* Compute the determinant of the Jacobian matrix */
        double detDF = Jrr * Jtt - Jrt * Jtr;
        value = 0.25 * (h1 + h2) * (k1 + k2) * fabs(detDF);
    }
    else if (i_r == 0 && DirBC_Interior)
    {
        value = 1.0;
    }
    else if (i_r == grid->nr() - 1)
    {
        value = 1.0;
    }

    rhs_f[grid->index(i_r, i_theta)] *= value;
}

void GMGPolar::discretize_rhs_f(const Level& level, GPU_Vector<double>& rhs_f)
{
    const PolarGrid& grid = level.grid();
    assert(rhs_f.size() == grid.numberOfNodes());

    const GPU_Vector<double>& sin_theta_cache = level.levelCache().GPU_sin_theta();
    const GPU_Vector<double>& cos_theta_cache = level.levelCache().GPU_cos_theta();

    dim3 threadsPerBlock(16, 16);
    dim3 numBlocks((grid.nr() + threadsPerBlock.x - 1) / threadsPerBlock.x,
                   (grid.ntheta() + threadsPerBlock.y - 1) / threadsPerBlock.y);

    DomainGeometry* device_domain_geometry;
    cudaMalloc(&device_domain_geometry, sizeof(DomainGeometry));
    cudaMemcpy(device_domain_geometry, domain_geometry_.get(), sizeof(DomainGeometry), cudaMemcpyHostToDevice);

    discretize_rhs_f_kernel<<<numBlocks, threadsPerBlock>>>(
        level.device_grid(), rhs_f.data(), sin_theta_cache.data(), cos_theta_cache.data(), DirBC_Interior_, device_domain_geometry);

    cudaDeviceSynchronize();

    cudaFree(device_domain_geometry);
}
