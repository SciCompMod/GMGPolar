#include "../../include/CoarseSolver/coarseSolver.h"


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



#define GIVE_BOUNDARY_SYMMETRY_SHIFT(i_theta) \
do { \
    theta = grid_.theta(i_theta); \
    sin_theta = sin_theta_[i_theta]; \
    cos_theta = cos_theta_[i_theta]; \
    int i_r; \
    double h1, h2, k1, k2; \
    double coeff1, coeff2; \
    if(DirBC_Interior_){ \
        /* -------------------------- */ \
        /* Node on the inner boundary */ \
        /* -------------------------- */ \
        i_r = 0; \
        r = grid_.radius(i_r); \
        coeff_alpha = system_parameters_.alpha(r); \
        coeff_beta = system_parameters_.beta(r); \
        ARR_ATT_ART(domain_geometry_, r, theta, sin_theta, cos_theta, coeff_alpha, arr, att, art, detDF); \
        h2 = grid_.r_dist(i_r); \
        k1 = grid_.theta_dist(i_theta-1); \
        k2 = grid_.theta_dist(i_theta); \
        coeff2 = 0.5*(k1+k2)/h2; \
        /* Fill x(i+1,j) */ \
        x[grid_.index(i_r+1,i_theta)] -= \
            - coeff2 * arr * x[grid_.index(i_r,i_theta)] /* Left */ \
            + 0.25 * art * x[grid_.index(i_r,i_theta+1)] /* Top Left */ \
            - 0.25 * art * x[grid_.index(i_r,i_theta-1)]; /* Bottom Left */ \
        /* --------------------------- */ \
        /* Node next to inner boundary */ \
        /* --------------------------- */ \
        i_r = 1; \
        r = grid_.radius(i_r); \
        coeff_alpha = system_parameters_.alpha(r); \
        coeff_beta = system_parameters_.beta(r); \
        ARR_ATT_ART(domain_geometry_, r, theta, sin_theta, cos_theta, coeff_alpha, \
                arr, att, art, detDF); \
        h1 = grid_.r_dist(i_r-1); \
        k1 = grid_.theta_dist(i_theta-1); \
        k2 = grid_.theta_dist(i_theta); \
        coeff1 = 0.5*(k1+k2)/h1; \
        /* Fill x(i,j) */ \
        x[grid_.index(i_r,i_theta)] -= \
            - coeff1 * arr * x[grid_.index(i_r-1,i_theta)]; /* Left */ \
        /* Fill x(i,j-1) */ \
        x[grid_.index(i_r,i_theta-1)] -= \
            + 0.25 * art * x[grid_.index(i_r-1,i_theta)]; /* Top Left */ \
        /* Fill x(i,j+1) */ \
        x[grid_.index(i_r,i_theta+1)] -= \
            - 0.25 * art * x[grid_.index(i_r-1,i_theta)]; /* Bottom Left */ \
    } \
    /* --------------------------- */ \
    /* Node next to outer boundary */ \
    /* --------------------------- */ \
    i_r = grid_.nr()-2; \
    r = grid_.radius(i_r); \
    coeff_alpha = system_parameters_.alpha(r); \
    coeff_beta = system_parameters_.beta(r); \
    ARR_ATT_ART(domain_geometry_, r, theta, sin_theta, cos_theta, coeff_alpha, \
            arr, att, art, detDF); \
    h2 = grid_.r_dist(i_r); \
    k1 = grid_.theta_dist(i_theta-1); \
    k2 = grid_.theta_dist(i_theta); \
    coeff2 = 0.5*(k1+k2)/h2; \
    /* Fill result(i,j) */ \
    x[grid_.index(i_r,i_theta)] -=  \
        - coeff2 * arr * x[grid_.index(i_r+1,i_theta)]; /* Right */ \
    /* Fill result(i,j-1) */ \
    x[grid_.index(i_r,i_theta-1)] -= \
        - 0.25 * art * x[grid_.index(i_r+1,i_theta)]; /* Top Right */ \
    /* Fill result(i,j+1) */ \
    x[grid_.index(i_r,i_theta+1)] -= \
        + 0.25 * art * x[grid_.index(i_r+1,i_theta)]; /* Bottom Right */ \
    /* -------------------------- */ \
    /* Node on the outer boundary */ \
    /* -------------------------- */ \
    i_r = grid_.nr()-1; \
    r = grid_.radius(i_r); \
    coeff_alpha = system_parameters_.alpha(r); \
    coeff_beta = system_parameters_.beta(r); \
    ARR_ATT_ART(domain_geometry_, r, theta, sin_theta, cos_theta, coeff_alpha, \
            arr, att, art, detDF); \
    h1 = grid_.r_dist(i_r-1); \
    k1 = grid_.theta_dist(i_theta-1); \
    k2 = grid_.theta_dist(i_theta); \
    coeff1 = 0.5*(k1+k2)/h1; \
    /* Fill result(i-1,j) */ \
    x[grid_.index(i_r-1,i_theta)] -= \
        - coeff1 * arr * x[grid_.index(i_r,i_theta)] /* Right */ \
        - 0.25 * art * x[grid_.index(i_r,i_theta+1)] /* Top Right */ \
        + 0.25 * art * x[grid_.index(i_r,i_theta-1)]; /* Bottom Right */ \
} while(0)


/* ----------------------- */
/* Boundary Symmetry Shift */
/* ----------------------- */

void CoarseSolver::subtractSymmetryShift(Vector<double>& x){
    assert(x.size() == grid_.number_of_nodes());

    omp_set_num_threads(maxOpenMPThreads_);

    const int additionalBoundaryTasks = grid_.ntheta() % 3;
    const int numBoundaryTasks = grid_.ntheta() - additionalBoundaryTasks;

    assert(numBoundaryTasks >= 3 && numBoundaryTasks % 3 == 0);

    /* Make sure to deallocate at the end */
    int* dep = new int[numBoundaryTasks];

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
            /* -------------- */
            /* Boundary Tasks */
            /* -------------- */

            /* Mod 0 Boundary */
            for(int boundary_task = 0; boundary_task < numBoundaryTasks; boundary_task += 3) {
                #pragma omp task \
                    depend(out: dep[boundary_task])
                {
                    if(boundary_task > 0){
                        int i_theta = boundary_task + additionalBoundaryTasks;    
                        GIVE_BOUNDARY_SYMMETRY_SHIFT(i_theta);
                    } else{
                        if(additionalBoundaryTasks == 0){
                            GIVE_BOUNDARY_SYMMETRY_SHIFT(0);
                        } 
                        else if(additionalBoundaryTasks >= 1){
                            GIVE_BOUNDARY_SYMMETRY_SHIFT(0);
                            GIVE_BOUNDARY_SYMMETRY_SHIFT(1);
                        }
                    }
                }
            }
            /* Mod 1 Boundary */
            for(int boundary_task = 1; boundary_task < numBoundaryTasks; boundary_task += 3) {
                #pragma omp task \
                    depend(out: dep[boundary_task]) \
                    depend(in: dep[boundary_task-1], dep[(boundary_task+2) % numBoundaryTasks])   
                {
                    if(boundary_task > 0){
                        int i_theta = boundary_task + additionalBoundaryTasks;    
                        GIVE_BOUNDARY_SYMMETRY_SHIFT(i_theta);
                    } else {
                        if(additionalBoundaryTasks == 0){
                            GIVE_BOUNDARY_SYMMETRY_SHIFT(1);
                        } 
                        else if(additionalBoundaryTasks == 1){
                            GIVE_BOUNDARY_SYMMETRY_SHIFT(2);
                        }
                        else if(additionalBoundaryTasks == 2){
                            GIVE_BOUNDARY_SYMMETRY_SHIFT(2);
                            GIVE_BOUNDARY_SYMMETRY_SHIFT(3);
                        }
                    }
                }
            }
            /* Mod 2 Boundary */
            for(int boundary_task = 2; boundary_task < numBoundaryTasks; boundary_task += 3) {
                #pragma omp task \
                    depend(out: dep[boundary_task]) \
                    depend(in: dep[boundary_task-1], dep[(boundary_task+2) % numBoundaryTasks])   
                {
                    int i_theta = boundary_task + additionalBoundaryTasks;    
                    GIVE_BOUNDARY_SYMMETRY_SHIFT(i_theta);
                }
            }
        }
    }
    omp_set_num_threads(maxOpenMPThreads_);

    delete[] dep;
}
