#include "../../include/Residual/residual.h"

Residual::Residual(const PolarGrid& grid, const LevelCache& level_cache, 
                    const DomainGeometry& domain_geometry,
                    bool DirBC_Interior, int num_omp_threads
) :
    grid_(grid),
    sin_theta_cache_(level_cache.sin_theta()),
    cos_theta_cache_(level_cache.cos_theta()),
    coeff_alpha_cache_(level_cache.coeff_alpha()),
    coeff_beta_cache_(level_cache.coeff_beta()),
    domain_geometry_(domain_geometry),
    DirBC_Interior_(DirBC_Interior),
    num_omp_threads_(num_omp_threads)
{}


#define COMPUTE_JACOBIAN_ELEMENTS(domain_geometry, r, theta, sin_theta, cos_theta, coeff_alpha, \
    arr, att, art, detDF) \
do { \
    /* Calculate the elements of the Jacobian matrix for the transformation mapping */ \
    /* The Jacobian matrix is: */ \
    /* [Jrr, Jrt] */ \
    /* [Jtr, Jtt] */ \
    const double Jrr = domain_geometry.dFx_dr(r, theta, sin_theta, cos_theta); \
    const double Jtr = domain_geometry.dFy_dr(r, theta, sin_theta, cos_theta); \
    const double Jrt = domain_geometry.dFx_dt(r, theta, sin_theta, cos_theta); \
    const double Jtt = domain_geometry.dFy_dt(r, theta, sin_theta, cos_theta); \
    /* Compute the determinant of the Jacobian matrix */ \
    detDF = Jrr * Jtt - Jrt * Jtr; \
    /* Compute the elements of the symmetric matrix: */ \
    /* 0.5 * alpha * DF^{-1} * DF^{-T} * |det(DF)| */ \
    /* which is represented by: */ \
    /* [arr, 0.5*art] */ \
    /* [0.5*atr, att] */ \
    arr = 0.5 * (Jtt * Jtt + Jrt * Jrt) * coeff_alpha / fabs(detDF); \
    att = 0.5 * (Jtr * Jtr + Jrr * Jrr) * coeff_alpha / fabs(detDF); \
    art = (- Jtt * Jtr - Jrt * Jrr) * coeff_alpha / fabs(detDF); \
    /* Note that the inverse Jacobian matrix DF^{-1} is: */ \
    /* 1.0 / det(DF) *   */ \
    /* [Jtt, -Jrt] */ \
    /* [-Jtr, Jrr] */ \
} while(0) \


#define NODE_APPLY_A_GIVE(i_r, i_theta, r, theta, sin_theta, cos_theta, \
    grid, DirBC_Interior, \
    result, x, factor, \
    arr, att, art, detDF, coeff_beta) \
do { \
    /* -------------------- */ \
    /* Node in the interior */ \
    /* -------------------- */ \
    if (i_r > 1 && i_r < grid.nr() - 2) { \
        double h1 = grid.radialSpacing(i_r-1); \
        double h2 = grid.radialSpacing(i_r); \
        double k1 = grid.angularSpacing(i_theta-1); \
        double k2 = grid.angularSpacing(i_theta); \
        double coeff1 = 0.5*(k1+k2)/h1; \
        double coeff2 = 0.5*(k1+k2)/h2; \
        double coeff3 = 0.5*(h1+h2)/k1; \
        double coeff4 = 0.5*(h1+h2)/k2; \
        /* Fill result(i,j) */ \
        result[grid.index(i_r,i_theta)] += factor * ( \
            0.25 * (h1+h2)*(k1+k2) * coeff_beta * fabs(detDF) * x[grid.index(i_r,i_theta)] /* beta_{i,j} */ \
            - coeff1 * arr * x[grid.index(i_r-1,i_theta)] /* Left */ \
            - coeff2 * arr * x[grid.index(i_r+1,i_theta)] /* Right */ \
            - coeff3 * att * x[grid.index(i_r,i_theta-1)] /* Bottom */ \
            - coeff4 * att * x[grid.index(i_r,i_theta+1)] /* Top */ \
            /* Center: (Left, Right, Bottom, Top) */ \
            + ((coeff1 + coeff2) * arr + (coeff3 + coeff4) * att) * x[grid.index(i_r,i_theta)] ); \
        /* Fill result(i-1,j) */ \
        result[grid.index(i_r-1,i_theta)] += factor * ( \
            - coeff1 * arr * x[grid.index(i_r,i_theta)] /* Right */ \
            + coeff1 * arr * x[grid.index(i_r-1,i_theta)] /* Center: (Right) */ \
            - 0.25 * art * x[grid.index(i_r,i_theta+1)] /* Top Right */ \
            + 0.25 * art * x[grid.index(i_r,i_theta-1)] ); /* Bottom Right */ \
        /* Fill result(i+1,j) */ \
        result[grid.index(i_r+1,i_theta)] += factor * ( \
            - coeff2 * arr * x[grid.index(i_r,i_theta)] /* Left */ \
            + coeff2 * arr * x[grid.index(i_r+1,i_theta)] /* Center: (Left) */ \
            + 0.25 * art * x[grid.index(i_r,i_theta+1)] /* Top Left */ \
            - 0.25 * art * x[grid.index(i_r,i_theta-1)] ); /* Bottom Left */ \
        /* Fill result(i,j-1) */ \
        result[grid.index(i_r,i_theta-1)] += factor * ( \
            - coeff3 * att * x[grid.index(i_r,i_theta)] /* Top */ \
            + coeff3 * att * x[grid.index(i_r,i_theta-1)] /* Center: (Top) */ \
            - 0.25 * art * x[grid.index(i_r+1,i_theta)] /* Top Right */ \
            + 0.25 * art * x[grid.index(i_r-1,i_theta)] ); /* Top Left */ \
        /* Fill result(i,j+1) */ \
        result[grid.index(i_r,i_theta+1)] += factor * ( \
            - coeff4 * att * x[grid.index(i_r,i_theta)] /* Bottom */ \
            + coeff4 * att * x[grid.index(i_r,i_theta+1)] /* Center: (Bottom) */ \
            + 0.25 * art * x[grid.index(i_r+1,i_theta)] /* Bottom Right */ \
            - 0.25 * art * x[grid.index(i_r-1,i_theta)] ); /* Bottom Left */ \
    /* -------------------------- */ \
    /* Node on the inner boundary */ \
    /* -------------------------- */ \
    } else if (i_r == 0) { \
        /* ------------------------------------------------ */ \
        /* Case 1: Dirichlet boundary on the inner boundary */ \
        /* ------------------------------------------------ */ \
        if(DirBC_Interior){ \
            /* Fill result(i,j) */ \
            result[grid.index(i_r,i_theta)] += factor * x[grid.index(i_r,i_theta)]; \
            /* Give value to the interior nodes! */ \
            double h2 = grid.radialSpacing(i_r); \
            double k1 = grid.angularSpacing(i_theta-1); \
            double k2 = grid.angularSpacing(i_theta); \
            double coeff2 = 0.5*(k1+k2)/h2; \
            /* Fill result(i+1,j) */ \
            result[grid.index(i_r+1,i_theta)] += factor * ( \
                - coeff2 * arr * x[grid.index(i_r,i_theta)] /* Left */ \
                + coeff2 * arr * x[grid.index(i_r+1,i_theta)] /* Center: (Left) */ \
                + 0.25 * art * x[grid.index(i_r,i_theta+1)] /* Top Left */ \
                - 0.25 * art * x[grid.index(i_r,i_theta-1)] ); /* Bottom Left */ \
        } else{ \
            /* ------------------------------------------------------------- */ \
            /* Case 2: Across origin discretization on the interior boundary */ \
            /* ------------------------------------------------------------- */ \
            /* h1 gets replaced with 2 * R0. */ \
            /* (i_r-1,i_theta) gets replaced with (i_r, i_theta + (grid.ntheta()>>1)). */ \
            /* Some more adjustments from the changing the 9-point stencil to the artifical 7-point stencil. */ \
            double h1 = 2.0 * grid.radius(0); \
            double h2 = grid.radialSpacing(i_r); \
            double k1 = grid.angularSpacing(i_theta-1); \
            double k2 = grid.angularSpacing(i_theta); \
            double coeff1 = 0.5*(k1+k2)/h1; \
            double coeff2 = 0.5*(k1+k2)/h2; \
            double coeff3 = 0.5*(h1+h2)/k1; \
            double coeff4 = 0.5*(h1+h2)/k2; \
            /* Fill result(i,j) */ \
            result[grid.index(i_r,i_theta)] += factor * ( \
                0.25 * (h1+h2)*(k1+k2) * coeff_beta * fabs(detDF) * x[grid.index(i_r,i_theta)] /* beta_{i,j} */ \
                - coeff1 * arr * x[grid.index(i_r, i_theta + (grid.ntheta()>>1))] /* Left */ \
                - coeff2 * arr * x[grid.index(i_r+1,i_theta)] /* Right */ \
                - coeff3 * att * x[grid.index(i_r,i_theta-1)] /* Bottom */ \
                - coeff4 * att * x[grid.index(i_r,i_theta+1)] /* Top */ \
                /* Center: (Left, Right, Bottom, Top) */ \
                + ((coeff1 + coeff2) * arr + (coeff3 + coeff4) * att) * x[grid.index(i_r,i_theta)] ); \
            /* Fill result(i-1,j) */ \
            /* From view the view of the across origin node, the directions are roatated by 180 degrees in the stencil! */ \
            result[grid.index(i_r, i_theta + (grid.ntheta()>>1))] += factor * ( \
                - coeff1 * arr * x[grid.index(i_r,i_theta)] /* Right -> Left */ \
                + coeff1 * arr * x[grid.index(i_r, i_theta + (grid.ntheta()>>1))] ); /* Center: (Right) -> Center: (Left)*/ \
            /*  + 0.25 * art * x[grid.index(i_r,i_theta+1)]; // Top Right -> Bottom Left: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */ \
            /*  - 0.25 * art * x[grid.index(i_r,i_theta-1)]; // Bottom Right -> Top Left: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */ \
            /* Fill result(i+1,j) */ \
            result[grid.index(i_r+1,i_theta)] += factor * ( \
                - coeff2 * arr * x[grid.index(i_r,i_theta)] /* Left */ \
                + coeff2 * arr * x[grid.index(i_r+1,i_theta)] /* Center: (Left) */ \
                + 0.25 * art * x[grid.index(i_r,i_theta+1)] /* Top Left */ \
                - 0.25 * art * x[grid.index(i_r,i_theta-1)] ); /* Bottom Left */ \
            /* Fill result(i,j-1) */ \
            result[grid.index(i_r,i_theta-1)] += factor * ( \
                - coeff3 * att * x[grid.index(i_r,i_theta)] /* Top */ \
                + coeff3 * att * x[grid.index(i_r,i_theta-1)] /* Center: (Top) */ \
                - 0.25 * art * x[grid.index(i_r+1,i_theta)] ); /* Top Right */ \
            /*  + 0.25 * art * x[grid.index(i_r-1,i_theta)]; // Top Left: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */ \
            /* Fill result(i,j+1) */ \
            result[grid.index(i_r,i_theta+1)] += factor * ( \
                - coeff4 * att * x[grid.index(i_r,i_theta)] /* Bottom */ \
                + coeff4 * att * x[grid.index(i_r,i_theta+1)] /* Center: (Bottom) */ \
                + 0.25 * art * x[grid.index(i_r+1,i_theta)] ); /* Bottom Right */ \
            /*  - 0.25 * art * x[grid.index(i_r-1,i_theta)]; // Bottom Left: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */ \
        } \
    /* ------------------------------- */ \
    /* Node next to the inner boundary */ \
    /* ------------------------------- */ \
    } else if (i_r == 1) { \
        double h1 = grid.radialSpacing(i_r-1); \
        double h2 = grid.radialSpacing(i_r); \
        double k1 = grid.angularSpacing(i_theta-1); \
        double k2 = grid.angularSpacing(i_theta); \
        double coeff1 = 0.5*(k1+k2)/h1; \
        double coeff2 = 0.5*(k1+k2)/h2; \
        double coeff3 = 0.5*(h1+h2)/k1; \
        double coeff4 = 0.5*(h1+h2)/k2; \
        /* Fill result(i,j) */ \
        result[grid.index(i_r,i_theta)] += factor * ( \
            0.25 * (h1+h2)*(k1+k2) * coeff_beta * fabs(detDF) * x[grid.index(i_r,i_theta)] /* beta_{i,j} */ \
            - coeff1 * arr * x[grid.index(i_r-1,i_theta)] /* Left */ \
            - coeff2 * arr * x[grid.index(i_r+1,i_theta)] /* Right */ \
            - coeff3 * att * x[grid.index(i_r,i_theta-1)] /* Bottom */ \
            - coeff4 * att * x[grid.index(i_r,i_theta+1)] /* Top */ \
            /* Center: (Left, Right, Bottom, Top) */ \
            + ((coeff1 + coeff2) * arr + (coeff3 + coeff4) * att) * x[grid.index(i_r,i_theta)] ); \
        /* Fill result(i-1,j) */ \
        if(!DirBC_Interior){ /* Don't give to the inner dirichlet boundary! */ \
            result[grid.index(i_r-1,i_theta)] += factor * ( \
                - coeff1 * arr * x[grid.index(i_r,i_theta)] /* Right */ \
                + coeff1 * arr * x[grid.index(i_r-1,i_theta)] /* Center: (Right) */ \
                - 0.25 * art * x[grid.index(i_r,i_theta+1)] /* Top Right */ \
                + 0.25 * art * x[grid.index(i_r,i_theta-1)] ); /* Bottom Right */ \
        } \
        /* Fill result(i+1,j) */ \
        result[grid.index(i_r+1,i_theta)] += factor * ( \
            - coeff2 * arr * x[grid.index(i_r,i_theta)] /* Left */ \
            + coeff2 * arr * x[grid.index(i_r+1,i_theta)] /* Center: (Left) */ \
            + 0.25 * art * x[grid.index(i_r,i_theta+1)] /* Top Left */ \
            - 0.25 * art * x[grid.index(i_r,i_theta-1)] ); /* Bottom Left */ \
        /* Fill result(i,j-1) */ \
        result[grid.index(i_r,i_theta-1)] += factor * ( \
            - coeff3 * att * x[grid.index(i_r,i_theta)] /* Top */ \
            + coeff3 * att * x[grid.index(i_r,i_theta-1)] /* Center: (Top) */ \
            - 0.25 * art * x[grid.index(i_r+1,i_theta)] /* Top Right */ \
            + 0.25 * art * x[grid.index(i_r-1,i_theta)] ); /* Top Left */ \
        /* Fill result(i,j+1) */ \
        result[grid.index(i_r,i_theta+1)] += factor * ( \
            - coeff4 * att * x[grid.index(i_r,i_theta)] /* Bottom */ \
            + coeff4 * att * x[grid.index(i_r,i_theta+1)] /* Center: (Bottom) */ \
            + 0.25 * art * x[grid.index(i_r+1,i_theta)] /* Bottom Right */ \
            - 0.25 * art * x[grid.index(i_r-1,i_theta)] ); /* Bottom Left */ \
    /* ------------------------------- */ \
    /* Node next to the outer boundary */ \
    /* ------------------------------- */ \
    } else if (i_r == grid.nr() - 2) { \
        double h1 = grid.radialSpacing(i_r-1); \
        double h2 = grid.radialSpacing(i_r); \
        double k1 = grid.angularSpacing(i_theta-1); \
        double k2 = grid.angularSpacing(i_theta); \
        double coeff1 = 0.5*(k1+k2)/h1; \
        double coeff2 = 0.5*(k1+k2)/h2; \
        double coeff3 = 0.5*(h1+h2)/k1; \
        double coeff4 = 0.5*(h1+h2)/k2; \
        /* Fill result(i,j) */ \
        result[grid.index(i_r,i_theta)] += factor * ( \
            0.25 * (h1+h2)*(k1+k2) * coeff_beta * fabs(detDF) * x[grid.index(i_r,i_theta)] /* beta_{i,j} */ \
            - coeff1 * arr * x[grid.index(i_r-1,i_theta)] /* Left */ \
            - coeff2 * arr * x[grid.index(i_r+1,i_theta)] /* Right */ \
            - coeff3 * att * x[grid.index(i_r,i_theta-1)] /* Bottom */ \
            - coeff4 * att * x[grid.index(i_r,i_theta+1)] /* Top */ \
            /* Center: (Left, Right, Bottom, Top) */ \
            + ((coeff1 + coeff2) * arr + (coeff3 + coeff4) * att) * x[grid.index(i_r,i_theta)] ); \
        /* Fill result(i-1,j) */ \
        result[grid.index(i_r-1,i_theta)] += factor * ( \
            - coeff1 * arr * x[grid.index(i_r,i_theta)] /* Right */ \
            + coeff1 * arr * x[grid.index(i_r-1,i_theta)] /* Center: (Right) */ \
            - 0.25 * art * x[grid.index(i_r,i_theta+1)] /* Top Right */ \
            + 0.25 * art * x[grid.index(i_r,i_theta-1)] ); /* Bottom Right */ \
        /* Don't give to the outer dirichlet boundary! */ \
        /* Fill result(i+1,j) */ \
        /* result[grid.index(i_r+1,i_theta)] += factor * ( */ \
        /*     - coeff2 * arr * x[grid.index(i_r,i_theta)] // Left */ \
        /*     + coeff2 * arr * x[grid.index(i_r+1,i_theta)] // Center: (Left) */ \
        /*     + 0.25 * art * x[grid.index(i_r,i_theta+1)] // Top Left */ \
        /*     - 0.25 * art * x[grid.index(i_r,i_theta-1)] ); // Bottom Left */ \
        /* Fill result(i,j-1) */ \
        result[grid.index(i_r,i_theta-1)] += factor * ( \
            - coeff3 * att * x[grid.index(i_r,i_theta)] /* Top */ \
            + coeff3 * att * x[grid.index(i_r,i_theta-1)] /* Center: (Top) */ \
            - 0.25 * art * x[grid.index(i_r+1,i_theta)] /* Top Right */ \
            + 0.25 * art * x[grid.index(i_r-1,i_theta)] ); /* Top Left */ \
        /* Fill result(i,j+1) */ \
        result[grid.index(i_r,i_theta+1)] += factor * ( \
            - coeff4 * att * x[grid.index(i_r,i_theta)] /* Bottom */ \
            + coeff4 * att * x[grid.index(i_r,i_theta+1)] /* Center: (Bottom) */ \
            + 0.25 * art * x[grid.index(i_r+1,i_theta)] /* Bottom Right */ \
            - 0.25 * art * x[grid.index(i_r-1,i_theta)] ); /* Bottom Left */ \
    /* ----------------------------- */ \
    /* Node on to the outer boundary */ \
    /* ----------------------------- */ \
    } else if (i_r == grid.nr() - 1) { \
        /* Fill result of (i,j) */ \
        result[grid.index(i_r,i_theta)] += factor * x[grid.index(i_r,i_theta)]; \
        /* Give value to the interior nodes! */ \
        double h1 = grid.radialSpacing(i_r-1); \
        double k1 = grid.angularSpacing(i_theta-1); \
        double k2 = grid.angularSpacing(i_theta); \
        double coeff1 = 0.5*(k1+k2)/h1; \
        /* Fill result(i-1,j) */ \
        result[grid.index(i_r-1,i_theta)] += factor * ( \
            - coeff1 * arr * x[grid.index(i_r,i_theta)] /* Right */ \
            + coeff1 * arr * x[grid.index(i_r-1,i_theta)] /* Center: (Right) */ \
            - 0.25 * art * x[grid.index(i_r,i_theta+1)] /* Top Right */ \
            + 0.25 * art * x[grid.index(i_r,i_theta-1)] ); /* Bottom Right */ \
    } \
} while(0)


void Residual::applyAGiveCircleSection(const int i_r, Vector<double>& result, const Vector<double>& x, const double& factor) const 
{
    const double r = grid_.radius(i_r);
    const double coeff_alpha = coeff_alpha_cache_[i_r];
    const double coeff_beta = coeff_beta_cache_[i_r];
    for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++){
        const double theta = grid_.theta(i_theta);
        const double sin_theta = sin_theta_cache_[i_theta];
        const double cos_theta = cos_theta_cache_[i_theta];
        /* Compute arr, att, art, detDF value at the current node */
        double arr, att, art, detDF;
        COMPUTE_JACOBIAN_ELEMENTS(domain_geometry_, r, theta, sin_theta, cos_theta, coeff_alpha, 
            arr, att, art, detDF);
        NODE_APPLY_A_GIVE(i_r, i_theta, r, theta, sin_theta, cos_theta,
            grid_, DirBC_Interior_,
            result, x, factor,
            arr, att, art, detDF, coeff_beta);
    }
}


void Residual::applyAGiveRadialSection(const int i_theta, Vector<double>& result, const Vector<double>& x, const double& factor) const 
{
    const double theta = grid_.theta(i_theta);
    const double sin_theta = sin_theta_cache_[i_theta];
    const double cos_theta = cos_theta_cache_[i_theta];
    for (int i_r = grid_.numberSmootherCircles(); i_r < grid_.nr(); i_r++){
        const double r = grid_.radius(i_r);
        const double coeff_alpha = coeff_alpha_cache_[i_r];
        const double coeff_beta = coeff_beta_cache_[i_r];
        // Compute arr, att, art, detDF value at the current node 
        double arr, att, art, detDF;
        COMPUTE_JACOBIAN_ELEMENTS(domain_geometry_, r, theta, sin_theta, cos_theta, coeff_alpha,
            arr, att, art, detDF);
        // Build solver matrix at the current node
        NODE_APPLY_A_GIVE(i_r, i_theta, r, theta, sin_theta, cos_theta,
            grid_, DirBC_Interior_,
            result, x, factor,
            arr, att, art, detDF, coeff_beta);
    }
}


/* ------------------ */
/* result = rhs - A*x */
void Residual::computeResidual(Vector<double>& result, const Vector<double>& rhs, const Vector<double>& x) const {
    assert(result.size() == x.size());

    bool use_simple_parallelism = true; // Fastest: true

    omp_set_num_threads(num_omp_threads_);

    result = rhs;

    const double factor = -1.0;

    if(omp_get_max_threads() == 1){
        /* Single-threaded execution */
        for(int i_r = 0; i_r < grid_.numberSmootherCircles(); i_r++) {
            applyAGiveCircleSection(i_r, result, x, factor);
        }
        for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++) {
            applyAGiveRadialSection(i_theta, result, x, factor);
        }
    }
    else{
        /* Multi-threaded execution */
        if(use_simple_parallelism) {
            /* For-Loop based parallelism */
            const int num_circle_tasks = grid_.numberSmootherCircles();
            const int additional_radial_tasks = grid_.ntheta() % 3;
            const int num_radial_tasks = grid_.ntheta() - additional_radial_tasks;

            #pragma omp parallel
            {
                #pragma omp for
                for(int circle_task = 0; circle_task < num_circle_tasks; circle_task += 3) {
                    int i_r = grid_.numberSmootherCircles() - circle_task - 1;   
                    applyAGiveCircleSection(i_r, result, x, factor);
                }
                #pragma omp for
                for(int circle_task = 1; circle_task < num_circle_tasks; circle_task += 3) {
                    int i_r = grid_.numberSmootherCircles() - circle_task - 1;   
                    applyAGiveCircleSection(i_r, result, x, factor);
                }
                #pragma omp for nowait
                for(int circle_task = 2; circle_task < num_circle_tasks; circle_task += 3) {
                    int i_r = grid_.numberSmootherCircles() - circle_task - 1;   
                    applyAGiveCircleSection(i_r, result, x, factor);
                }

                #pragma omp for
                for (int radial_task = 0; radial_task < num_radial_tasks; radial_task += 3) {
                    if(radial_task > 0){
                        int i_theta = radial_task + additional_radial_tasks;    
                        applyAGiveRadialSection(i_theta, result, x, factor);
                    } else{
                        if(additional_radial_tasks == 0){
                            applyAGiveRadialSection(0, result, x, factor);
                        } 
                        else if(additional_radial_tasks >= 1){
                            applyAGiveRadialSection(0, result, x, factor);
                            applyAGiveRadialSection(1, result, x, factor);
                        }
                    }
                }
                #pragma omp for
                for (int radial_task = 1; radial_task < num_radial_tasks; radial_task += 3) {
                    if(radial_task > 1){
                        int i_theta = radial_task + additional_radial_tasks;    
                        applyAGiveRadialSection(i_theta, result, x, factor);
                    } else {
                        if(additional_radial_tasks == 0){
                            applyAGiveRadialSection(1, result, x, factor);
                        } 
                        else if(additional_radial_tasks == 1){
                            applyAGiveRadialSection(2, result, x, factor);
                        }
                        else if(additional_radial_tasks == 2){
                            applyAGiveRadialSection(2, result, x, factor);
                            applyAGiveRadialSection(3, result, x, factor);
                        }
                    }
                }
                #pragma omp for
                for (int radial_task = 2; radial_task < num_radial_tasks; radial_task += 3) {
                    int i_theta = radial_task + additional_radial_tasks;    
                    applyAGiveRadialSection(i_theta, result, x, factor);
                }
            }
        }
        else{
            /* Task dependency based parallelism */
            const int num_circle_tasks = grid_.numberSmootherCircles();
            const int additional_radial_tasks = grid_.ntheta() % 3;
            const int num_radial_tasks = grid_.ntheta() - additional_radial_tasks;

            assert(num_circle_tasks >= 2);
            assert(num_radial_tasks >= 3 && num_radial_tasks % 3 == 0);

            /* Make sure to deallocate at the end */
            const int boundary_margin = 2; // Additional space to ensure safe access
            int* circle_dep = new int[num_circle_tasks + boundary_margin];
            int* radial_dep = new int[num_radial_tasks];

            #pragma omp parallel 
            {
                #pragma omp single
                {
                    /* ------------ */
                    /* Circle Tasks */
                    /* ------------ */

                    /* Mod 0 Circles */
                    for(int circle_task = 0; circle_task < num_circle_tasks; circle_task += 3) {
                        #pragma omp task \
                            depend(out: circle_dep[circle_task])
                        {
                            int i_r = grid_.numberSmootherCircles() - circle_task - 1;    
                            applyAGiveCircleSection(i_r, result, x, factor);
                        }
                    }
                    /* Mod 2 Circles */
                    for(int circle_task = 1; circle_task < num_circle_tasks; circle_task += 3) {
                        #pragma omp task \
                            depend(out: circle_dep[circle_task]) \
                            depend(in: circle_dep[circle_task-1], circle_dep[circle_task+2])   
                        {
                            int i_r = grid_.numberSmootherCircles() - circle_task - 1;    
                            applyAGiveCircleSection(i_r, result, x, factor);
                        }
                    }
                    /* Mod 2 Circles */
                    for(int circle_task = 2; circle_task < num_circle_tasks; circle_task += 3) {
                        #pragma omp task \
                            depend(out: circle_dep[circle_task]) \
                            depend(in: circle_dep[circle_task-1], circle_dep[circle_task+2])   
                        {
                            int i_r = grid_.numberSmootherCircles() - circle_task - 1;    
                            applyAGiveCircleSection(i_r, result, x, factor);
                        }
                    }

                    /* ------------ */
                    /* Radial Tasks */
                    /* ------------ */

                    /* Mod 0 Radials */
                    for(int radial_task = 0; radial_task < num_radial_tasks; radial_task += 3) {
                        #pragma omp task \
                            depend(out: radial_dep[radial_task]) \
                            depend(in: circle_dep[1]) /* Wait for Circle Smoother */
                        {
                            if(radial_task > 0){
                                int i_theta = radial_task + additional_radial_tasks;    
                                applyAGiveRadialSection(i_theta, result, x, factor);
                            } else{
                                if(additional_radial_tasks == 0){
                                    applyAGiveRadialSection(0, result, x, factor);
                                } 
                                else if(additional_radial_tasks >= 1){
                                    applyAGiveRadialSection(0, result, x, factor);
                                    applyAGiveRadialSection(1, result, x, factor);
                                }
                            }
                        }
                    }
                    /* Mod 1 Radials */
                    for(int radial_task = 1; radial_task < num_radial_tasks; radial_task += 3) {
                        #pragma omp task \
                            depend(out: radial_dep[radial_task]) \
                            depend(in: \
                                radial_dep[radial_task-1], \
                                radial_dep[(radial_task+2) % num_radial_tasks])   
                        {
                            if(radial_task > 1){
                                int i_theta = radial_task + additional_radial_tasks;    
                                applyAGiveRadialSection(i_theta, result, x, factor);
                            } else {
                                if(additional_radial_tasks == 0){
                                    applyAGiveRadialSection(1, result, x, factor);
                                } 
                                else if(additional_radial_tasks == 1){
                                    applyAGiveRadialSection(2, result, x, factor);
                                }
                                else if(additional_radial_tasks == 2){
                                    applyAGiveRadialSection(2, result, x, factor);
                                    applyAGiveRadialSection(3, result, x, factor);
                                }
                            }
                        }
                    }
                    /* Mod 2 Radials */
                    for(int radial_task = 2; radial_task < num_radial_tasks; radial_task += 3) {
                        #pragma omp task \
                            depend(in: \
                                radial_dep[radial_task-1], \
                                radial_dep[(radial_task+2) % num_radial_tasks])   
                        {
                            int i_theta = radial_task + additional_radial_tasks;    
                            applyAGiveRadialSection(i_theta, result, x, factor);
                        }
                    }
                }
            }
            delete[] circle_dep;
            delete[] radial_dep;
        }
    }
}