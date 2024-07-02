#include "../../include/Residual/residual.h"


Residual::Residual(const PolarGrid& grid, const LevelCache& level_data, 
    const DomainGeometry& domain_geometry, const SystemParameters& system_parameters, const bool DirBC_Interior, 
    const int maxOpenMPThreads, const int openMPTaskThreads
) :
    grid_(grid), 
    sin_theta_(level_data.sin_theta()),
    cos_theta_(level_data.cos_theta()),
    domain_geometry_(domain_geometry),
    system_parameters_(system_parameters),
    DirBC_Interior_(DirBC_Interior),
    maxOpenMPThreads_(maxOpenMPThreads),
    openMPTaskThreads_(openMPTaskThreads)
{}



#define ARR_ATT_ART(domain_geometry, r, theta, sin_theta, cos_theta, coeff_alpha, \
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



#define NODE_APPLY_RESIDUAL_GIVE(i_r, i_theta, r, theta, sin_theta, cos_theta, \
    system_parameters, grid, DirBC_Interior, \
    result, rhs, x, factor, \
    arr, att, art, coeff_beta, detDF) \
do { \
    /* -------------------- */ \
    /* Node in the interior */ \
    /* -------------------- */ \
    if (i_r > 1 && i_r < grid.nr() - 2) { \
        double h1 = grid.r_dist(i_r-1); \
        double h2 = grid.r_dist(i_r); \
        double k1 = grid.theta_dist(i_theta-1); \
        double k2 = grid.theta_dist(i_theta); \
        double coeff1 = 0.5*(k1+k2)/h1; \
        double coeff2 = 0.5*(k1+k2)/h2; \
        double coeff3 = 0.5*(h1+h2)/k1; \
        double coeff4 = 0.5*(h1+h2)/k2; \
        /* Fill result of (i,j) */ \
        result[grid.index(i_r,i_theta)] += rhs[grid.index(i_r,i_theta)]; \
        /* result[grid.index(i_r,i_theta)] += rhs[grid.index(i_r,i_theta)]; */ \
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
            /* Fill result of (i,j) */ \
            result[grid.index(i_r,i_theta)] += rhs[grid.index(i_r,i_theta)]; /* Contains u_D_Interior */ \
            /* result[grid.index(i_r,i_theta)] += rhs[grid.index(i_r,i_theta)]; // Contains u_D_Interior */ \
            /* Fill result(i,j) */ \
            result[grid.index(i_r,i_theta)] += factor * x[grid.index(i_r,i_theta)]; \
            /* Give value to the interior nodes! */ \
            double h2 = grid.r_dist(i_r); \
            double k1 = grid.theta_dist(i_theta-1); \
            double k2 = grid.theta_dist(i_theta); \
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
            double h2 = grid.r_dist(i_r); \
            double k1 = grid.theta_dist(i_theta-1); \
            double k2 = grid.theta_dist(i_theta); \
            double coeff1 = 0.5*(k1+k2)/h1; \
            double coeff2 = 0.5*(k1+k2)/h2; \
            double coeff3 = 0.5*(h1+h2)/k1; \
            double coeff4 = 0.5*(h1+h2)/k2; \
            /* Fill result of (i,j) */ \
            result[grid.index(i_r,i_theta)] += rhs[grid.index(i_r,i_theta)]; \
            /* result[grid.index(i_r,i_theta)] += rhs[grid.index(i_r,i_theta)]; */ \
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
        double h1 = grid.r_dist(i_r-1); \
        double h2 = grid.r_dist(i_r); \
        double k1 = grid.theta_dist(i_theta-1); \
        double k2 = grid.theta_dist(i_theta); \
        double coeff1 = 0.5*(k1+k2)/h1; \
        double coeff2 = 0.5*(k1+k2)/h2; \
        double coeff3 = 0.5*(h1+h2)/k1; \
        double coeff4 = 0.5*(h1+h2)/k2; \
        /* Fill result of (i,j) */ \
        result[grid.index(i_r,i_theta)] += rhs[grid.index(i_r,i_theta)]; \
        /* result[grid.index(i_r,i_theta)] += rhs[grid.index(i_r,i_theta)]; */ \
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
        double h1 = grid.r_dist(i_r-1); \
        double h2 = grid.r_dist(i_r); \
        double k1 = grid.theta_dist(i_theta-1); \
        double k2 = grid.theta_dist(i_theta); \
        double coeff1 = 0.5*(k1+k2)/h1; \
        double coeff2 = 0.5*(k1+k2)/h2; \
        double coeff3 = 0.5*(h1+h2)/k1; \
        double coeff4 = 0.5*(h1+h2)/k2; \
        /* Fill result of (i,j) */ \
        result[grid.index(i_r,i_theta)] += rhs[grid.index(i_r,i_theta)]; \
        /* result[grid.index(i_r,i_theta)] += rhs[grid.index(i_r,i_theta)]; */ \
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
        result[grid.index(i_r,i_theta)] += rhs[grid.index(i_r,i_theta)]; /* Contains u_D */ \
        /* result[grid.index(i_r,i_theta)] += rhs[grid.index(i_r,i_theta)]; // Contains u_D */ \
        /* Dirichlet boundary */ \
        result[grid.index(i_r,i_theta)] += factor * x[grid.index(i_r,i_theta)]; \
        /* Give value to the interior nodes! */ \
        double h1 = grid.r_dist(i_r-1); \
        double k1 = grid.theta_dist(i_theta-1); \
        double k2 = grid.theta_dist(i_theta); \
        double coeff1 = 0.5*(k1+k2)/h1; \
        /* Fill result(i-1,j) */ \
        result[grid.index(i_r-1,i_theta)] += factor * ( \
            - coeff1 * arr * x[grid.index(i_r,i_theta)] /* Right */ \
            + coeff1 * arr * x[grid.index(i_r-1,i_theta)] /* Center: (Right) */ \
            - 0.25 * art * x[grid.index(i_r,i_theta+1)] /* Top Right */ \
            + 0.25 * art * x[grid.index(i_r,i_theta-1)] ); /* Bottom Right */ \
    } \
} while(0)



#define CIRCLE_SECTION_APPLY_RESIDUAL_GIVE(i_r) \
do { \
    r = grid_.radius(i_r); \
    coeff_alpha = system_parameters_.alpha(r); \
    coeff_beta = system_parameters_.beta(r); \
    for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++){ \
        theta = grid_.theta(i_theta); \
        sin_theta = sin_theta_[i_theta]; \
        cos_theta = cos_theta_[i_theta]; \
        \
        ARR_ATT_ART(domain_geometry_, r, theta, sin_theta, cos_theta, coeff_alpha, \
            arr, att, art, detDF); \
        \
        NODE_APPLY_RESIDUAL_GIVE(i_r, i_theta, r, theta, sin_theta, cos_theta, \
            system_parameters_, grid_, DirBC_Interior_, \
            result, rhs, x, factor, \
            arr, att, art, coeff_beta, detDF); \
    } \
} while(0)



#define RADIAL_SECTION_APPLY_RESIDUAL_GIVE(i_theta) \
do { \
    theta = grid_.theta(i_theta); \
    sin_theta = sin_theta_[i_theta]; \
    cos_theta = cos_theta_[i_theta]; \
    for (int i_r = grid_.numberSmootherCircles(); i_r < grid_.nr(); i_r++){ \
        r = grid_.radius(i_r); \
        coeff_alpha = system_parameters_.alpha(r); \
        coeff_beta = system_parameters_.beta(r); \
        /* Get arr, att, art, detDF value at the current node */ \
        ARR_ATT_ART(domain_geometry_, r, theta, sin_theta, cos_theta, coeff_alpha, \
            arr, att, art, detDF); \
        \
        NODE_APPLY_RESIDUAL_GIVE(i_r, i_theta, r, theta, sin_theta, cos_theta, \
            system_parameters_, grid_, DirBC_Interior_, \
            result, rhs, x, factor, \
            arr, att, art, coeff_beta, detDF); \
    } \
} while(0)



/* ------------------------- */
/* result = rhs + factor * A*x */
void Residual::computeResidual_V1(Vector<double>& result, const Vector<double>& rhs, const Vector<double>& x) const {
    assert(result.size() == x.size());

    const double factor = -1.0;

    omp_set_num_threads(maxOpenMPThreads_);
    assign(result, 0.0);

    const int numCircleTasks = grid_.numberSmootherCircles();
    const int additionalRadialTasks = grid_.ntheta() % 3;
    const int numRadialTasks = grid_.ntheta() - additionalRadialTasks;

    assert(numCircleTasks >= 2);
    assert(numRadialTasks >= 3 && numRadialTasks % 3 == 0);

    /* Make sure to deallocate at the end */
    int* dep = new int[numCircleTasks + numRadialTasks];

    omp_set_num_threads(openMPTaskThreads_);
    #pragma omp parallel num_threads(openMPTaskThreads_) /* Outside variable are shared by default */
    {
        /* Define thread-local variables */
        double r, theta;
        double sin_theta, cos_theta;
        double arr, att, art;
        double coeff_alpha, coeff_beta;
        double detDF;

        #pragma omp single
        {
            /* ------------ */
            /* Circle Tasks */
            /* ------------ */

            /* Mod 0 Circles */
            for(int circle_task = 0; circle_task < numCircleTasks; circle_task += 3) {
                #pragma omp task \
                    depend(out: dep[circle_task])
                {
                    int i_r = grid_.numberSmootherCircles() - circle_task - 1;    
                    CIRCLE_SECTION_APPLY_RESIDUAL_GIVE(i_r);
                }
            }
            /* Mod 2 Circles */
            for(int circle_task = 1; circle_task < numCircleTasks; circle_task += 3) {
                #pragma omp task \
                    depend(out: dep[circle_task]) \
                    depend(in: dep[circle_task-1], dep[circle_task+2])   
                {
                    int i_r = grid_.numberSmootherCircles() - circle_task - 1;    
                    CIRCLE_SECTION_APPLY_RESIDUAL_GIVE(i_r);
                }
            }
            /* Mod 2 Circles */
            for(int circle_task = 2; circle_task < numCircleTasks; circle_task += 3) {
                #pragma omp task \
                    depend(out: dep[circle_task]) \
                    depend(in: dep[circle_task-1], dep[circle_task+2])   
                {
                    int i_r = grid_.numberSmootherCircles() - circle_task - 1;    
                    CIRCLE_SECTION_APPLY_RESIDUAL_GIVE(i_r);
                }
            }

            /* ------------ */
            /* Radial Tasks */
            /* ------------ */

            /* Mod 0 Radials */
            for(int radial_task = 0; radial_task < numRadialTasks; radial_task += 3) {
                #pragma omp task \
                    depend(out: dep[numCircleTasks+radial_task]) \
                    depend(in: dep[1]) /* Wait for Circle Smoother */
                {
                    if(radial_task > 0){
                        int i_theta = radial_task + additionalRadialTasks;    
                        RADIAL_SECTION_APPLY_RESIDUAL_GIVE(i_theta);
                    } else{
                        if(additionalRadialTasks == 0){
                            RADIAL_SECTION_APPLY_RESIDUAL_GIVE(0);
                        } 
                        else if(additionalRadialTasks >= 1){
                            RADIAL_SECTION_APPLY_RESIDUAL_GIVE(0);
                            RADIAL_SECTION_APPLY_RESIDUAL_GIVE(1);
                        }
                    }
                }
            }
            /* Mod 1 Radials */
            for(int radial_task = 1; radial_task < numRadialTasks; radial_task += 3) {
                #pragma omp task \
                    depend(out: dep[numCircleTasks + radial_task]) \
                    depend(in: \
                        dep[1], /* Wait for Circle Smoother */ \
                        dep[numCircleTasks + radial_task-1], \
                        dep[numCircleTasks + (radial_task+2) % numRadialTasks])   
                {
                    if(radial_task > 1){
                        int i_theta = radial_task + additionalRadialTasks;    
                        RADIAL_SECTION_APPLY_RESIDUAL_GIVE(i_theta);
                    } else {
                        if(additionalRadialTasks == 0){
                            RADIAL_SECTION_APPLY_RESIDUAL_GIVE(1);
                        } 
                        else if(additionalRadialTasks == 1){
                            RADIAL_SECTION_APPLY_RESIDUAL_GIVE(2);
                        }
                        else if(additionalRadialTasks == 2){
                            RADIAL_SECTION_APPLY_RESIDUAL_GIVE(2);
                            RADIAL_SECTION_APPLY_RESIDUAL_GIVE(3);
                        }
                    }
                }
            }
            /* Mod 2 Radials */
            for(int radial_task = 2; radial_task < numRadialTasks; radial_task += 3) {
                #pragma omp task \
                    depend(out: dep[numCircleTasks + radial_task]) \
                    depend(in: \
                        dep[1], /* Wait for Circle Smoother */ \
                        dep[numCircleTasks + radial_task-1], \
                        dep[numCircleTasks + (radial_task+2) % numRadialTasks])   
                {
                    int i_theta = radial_task + additionalRadialTasks;    
                    RADIAL_SECTION_APPLY_RESIDUAL_GIVE(i_theta);
                }
            }
        }
    }

    delete[] dep;
    omp_set_num_threads(maxOpenMPThreads_);
}



/* ------------------------- */
/* result = f + factor * A*x */
void Residual::computeResidual_V2(Vector<double>& result, const Vector<double>& rhs, const Vector<double>& x) const{
    assert(result.size() == x.size());

    const double factor = -1.0;

    omp_set_num_threads(maxOpenMPThreads_);
    assign(result, 0.0);

    const int numCircleTasks = grid_.numberSmootherCircles();
    const int additionalRadialTasks = grid_.ntheta() % 3;
    const int numRadialTasks = grid_.ntheta() - additionalRadialTasks;

    assert(numCircleTasks >= 2);
    assert(numRadialTasks >= 3 && numRadialTasks % 3 == 0);

    #pragma omp parallel num_threads(maxOpenMPThreads_) /* Outside variable are shared by default */
    {
        /* Define thread-local variables */
        double r, theta;
        double sin_theta, cos_theta;
        double arr, att, art;
        double coeff_alpha, coeff_beta;
        double detDF;

        /* ------------ */
        /* Circle Tasks */
        /* ------------ */

        /* Mod 0 Circles */
        #pragma omp for
        for (int circle_task = 0; circle_task < numCircleTasks; circle_task += 3){
            int i_r = grid_.numberSmootherCircles() - circle_task - 1;    
            CIRCLE_SECTION_APPLY_RESIDUAL_GIVE(i_r);
        }
        #pragma omp barrier
        /* Mod 1 Circles */
        #pragma omp for
        for (int circle_task = 1; circle_task < numCircleTasks; circle_task += 3){
            int i_r = grid_.numberSmootherCircles() - circle_task - 1;    
            CIRCLE_SECTION_APPLY_RESIDUAL_GIVE(i_r);
        }
        /* Mod 2 Circles */
        #pragma omp barrier
        #pragma omp for nowait
        for (int circle_task = 2; circle_task < numCircleTasks; circle_task += 3){
            int i_r = grid_.numberSmootherCircles() - circle_task - 1;    
            CIRCLE_SECTION_APPLY_RESIDUAL_GIVE(i_r);
        }

        /* ------------ */
        /* Radial Tasks */
        /* ------------ */
        
        /* Mod 0 Radials */
        #pragma omp for
        for (int radial_task = 0; radial_task < numRadialTasks; radial_task += 3){
            if(radial_task > 0){
                int i_theta = radial_task + additionalRadialTasks;    
                RADIAL_SECTION_APPLY_RESIDUAL_GIVE(i_theta);
            } else{
                if(additionalRadialTasks == 0){
                    RADIAL_SECTION_APPLY_RESIDUAL_GIVE(0);
                } 
                else if(additionalRadialTasks >= 1){
                    RADIAL_SECTION_APPLY_RESIDUAL_GIVE(0);
                    RADIAL_SECTION_APPLY_RESIDUAL_GIVE(1);
                }
            }
        }
        #pragma omp barrier
        /* Mod 1 Radials */
        #pragma omp for
        for (int radial_task = 1; radial_task < numRadialTasks; radial_task += 3){
            if(radial_task > 1){
                int i_theta = radial_task + additionalRadialTasks;    
                RADIAL_SECTION_APPLY_RESIDUAL_GIVE(i_theta);
            } else {
                if(additionalRadialTasks == 0){
                    RADIAL_SECTION_APPLY_RESIDUAL_GIVE(1);
                } 
                else if(additionalRadialTasks == 1){
                    RADIAL_SECTION_APPLY_RESIDUAL_GIVE(2);
                }
                else if(additionalRadialTasks == 2){
                    RADIAL_SECTION_APPLY_RESIDUAL_GIVE(2);
                    RADIAL_SECTION_APPLY_RESIDUAL_GIVE(3);
                }
            }
        }
        #pragma omp barrier
        /* Mod 2 Radials */
        #pragma omp for
        for (int radial_task = 2; radial_task < numRadialTasks; radial_task += 3){
            int i_theta = radial_task + additionalRadialTasks;    
            RADIAL_SECTION_APPLY_RESIDUAL_GIVE(i_theta);
        }
    }
}



/* ------------------------- */
/* result = f + factor * A*x */
void Residual::computeResidual_V3(Vector<double>& result, const Vector<double>& rhs, const Vector<double>& x) const{
    assert(result.size() == x.size());

    const double factor = -1.0;

    omp_set_num_threads(maxOpenMPThreads_);
    assign(result, 0.0);

    // ------------------------ //
    // Custom Task Distribution //
    // ------------------------ //
    const int numThreads = omp_get_max_threads();
    omp_set_num_threads(numThreads);
    const int minimalChunkSize = 4;
    const int zone = 2;

    // Distribute Tasks to each thread
    TaskDistribution CircleSmootherTasks(grid_.numberSmootherCircles(), minimalChunkSize, numThreads);
    TaskDistribution RadialSmootherTasks(grid_.ntheta(), minimalChunkSize, numThreads);

    #pragma omp parallel num_threads(maxOpenMPThreads_) /* Outside variable are shared by default */
    {
        /* Define thread-local variables */
        double r, theta;
        double sin_theta, cos_theta;
        double arr, att, art;
        double coeff_alpha, coeff_beta;
        double detDF;

        const int threadID = omp_get_thread_num();
        // ----------------------------------------------------------- //
        // Take care of the separation strips of the circular smoother //
        // ----------------------------------------------------------- //
        const int circle_task_start = CircleSmootherTasks.getStart(threadID);
        const int circle_task_end = CircleSmootherTasks.getEnd(threadID);
        const int circle_task_separation = std::min(circle_task_end - circle_task_start, zone);

        for (int circle_task = circle_task_start; circle_task < circle_task_start + circle_task_separation; circle_task++){
            int i_r = grid_.numberSmootherCircles() - circle_task - 1;    
            CIRCLE_SECTION_APPLY_RESIDUAL_GIVE(i_r);
        }
        
        #pragma omp barrier

        // -------------------------------------------------------- //
        // Take care of the separation strips of the radial smoother //
        // -------------------------------------------------------- //
        const int radial_task_start = RadialSmootherTasks.getStart(threadID);
        const int radial_task_end = RadialSmootherTasks.getEnd(threadID);
        const int radial_task_separation = std::min(radial_task_end-radial_task_start, zone);

        for (int radial_task = radial_task_start; radial_task < radial_task_start + radial_task_separation; radial_task++){
            int i_theta = radial_task;  
            RADIAL_SECTION_APPLY_RESIDUAL_GIVE(i_theta);
        }

        #pragma omp barrier

        // ------------------------------------------ //
        // Take care of the circular smoother section //
        // ------------------------------------------ //
        for (int circle_task = circle_task_start + circle_task_separation; circle_task < circle_task_end; circle_task++){
            int i_r = grid_.numberSmootherCircles() - circle_task - 1;    
            CIRCLE_SECTION_APPLY_RESIDUAL_GIVE(i_r);
        }

        // ---------------------------------------- //
        // Take care of the radial smoother section //
        // ---------------------------------------- //
        for (int radial_task = radial_task_start + radial_task_separation; radial_task < radial_task_end; radial_task++){
            int i_theta = radial_task;  
            RADIAL_SECTION_APPLY_RESIDUAL_GIVE(i_theta);
        }
    }
}