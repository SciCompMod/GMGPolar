#include <iostream>
#include <vector>
#include <cmath>
#include <thread>
#include <omp.h>
#include <mutex>

#include "../include/InputFunctions/domain_geometry.h"
#include "../include/InputFunctions/system_parameters.h"
#include "../include/InputFunctions/exact_solution.h"
#include "../include/common/constants.h"

#include "../include/GMGPolar/gmgpolar.h"

// #include "../include/TaskDistribution/TaskDistribution.h"


#include "../include/linear_algebra/vector.h"
#include "../include/linear_algebra/matrix.h"
#include "../include/linear_algebra/operations.h"
#include <chrono>

#include <iostream>
#include <vector>
#include <omp.h>
#include <thread>
#include <chrono>
#include <random>
#include <memory>






#include <stdio.h>


  

#include "mpi.h" 
#include "dmumps_c.h"   
#include <likwid.h>





// int main(int argc, char** argv) {

//     LIKWID_MARKER_INIT;

//     LIKWID_MARKER_START("create_grid_polar");

    // Perform matrix multiplication
    
    // JOB (integer) must be initialized by the user on all processors before a call to MUMPS. 
    // It controls the main actions taken by MUMPS. It is not altered by MUMPS. 
    // Possible values of JOB are:
    // • JOB= –1 initializes an instance of the package.
    // • JOB= –2 terminates an instance of the package.
    // • JOB= –3 save / restore feature: removes data saved to disk.
    // • JOB= –4 after factorization or solve phases, frees all MUMPS internal data structures except the ones from analysis.
    // • JOB= 1 performs the analysis phase.
    // • JOB= 2 performs the factorization phase.
    // • JOB= 3 computes the solution.
    // • JOB= 4 combines the actions of JOB= 1 with those of JOB= 2.
    // • JOB= 5 combines the actions of JOB= 2 and JOB= 3.
    // • JOB= 6 combines the actions of calls with JOB= 1, JOB= 2, and JOB= 3.
    // • JOB= 7 save / restore feature: saves MUMPS internal data to disk.
    // • JOB= 8 save / restore feature: restores MUMPS internal data from disk.
    // • JOB= 9 computes before the solution phase a possible distribution for the right-hand sides.


//     // MUMPS_INT n = 3; // Size of the matrix
//     // MUMPS_INT nz = 4; // Number of non-zero elements
//     // MUMPS_INT irn[4] = {1, 2, 3, 1}; // Row indices of non-zero elements
//     // MUMPS_INT jcn[4] = {1, 2, 3, 3}; // Column indices of non-zero elements
//     // double a[4] = {1.0, 1.0, 1.0, 6.0}; // Values of non-zero elements
//     // MUMPS_INT nz_rhs = 3; // Number of right-hand sides
//     // double rhs[3] = {1.0, 2.0, 3.0}; // Right-hand side vector
//     // MUMPS_INT lrhs = 3; // Leading dimension of rhs (should be >= nz_rhs)
//     // MUMPS_INT ierr;

//     // DMUMPS_STRUC_C mumps;
//     // mumps.job = -1; // Initialize MUMPS
//     // mumps.par = 1; // Sequential execution
//     // mumps.sym = 1; // General matrix
//     // mumps.comm_fortran = -987654; // Dummy MPI communicator

//     // dmumps_c(&mumps);

//     // mumps.job = 4; // Solve Ax = b
//     // mumps.n = n;
//     // mumps.nz = nz;
//     // mumps.irn = irn;
//     // mumps.jcn = jcn;
//     // mumps.a = a;
//     // mumps.nz_rhs = nz_rhs;
//     // mumps.rhs = rhs;
//     // mumps.lrhs = lrhs;

//     // dmumps_c(&mumps);

//     // mumps.job = 3;
//     // dmumps_c(&mumps);

//     // if (mumps.info[0] != 0) {
//     //     std::cerr << "Error solving the system: " << mumps.info[0] << std::endl;
//     //     return 1;
//     // }

//     // std::cout << "Solution:" << std::endl;
//     // for (int i = 0; i < n; ++i) {
//     //     std::cout << mumps.rhs[i] << std::endl;
//     // }





//     int n = 3;
//     std::vector<std::tuple<int, int, double>> entries = {
//         {1, 1, 1.0}, {2, 2, 1.0}, {3, 3, 1.0}, {1, 3, 6.0}
//     };
//     SparseMatrix matrix(n, n, entries);
//     matrix.is_symmetric(true);

//     DMUMPS_STRUC_C mumps;
//     mumps.job = -1;
//     mumps.par = 1;
//     mumps.sym = (matrix.is_symmetric() ? 1 : 0);
//     mumps.comm_fortran = -987654;
//     dmumps_c(&mumps);

//     mumps.job = 4; 
//     assert(matrix.rows() == matrix.cols());
//     mumps.n = matrix.rows();
//     mumps.nz = matrix.non_zero_size();
//     mumps.irn = matrix.row_indices_data();
//     mumps.jcn = matrix.column_indices_data();
//     mumps.a = matrix.values_data();
//     dmumps_c(&mumps);

//     int nz_rhs = 7; 
//     int nrhs = 2;
//     double rhs[6] = {1.0, 2.0, 3.0, 1.0, 2.0, 3.0}; 
//     int lrhs = 3; 
    
//     mumps.job = 3; 
//     mumps.nrhs = nrhs;
//     mumps.nz_rhs = nz_rhs;
//     mumps.rhs = rhs;
//     mumps.lrhs = lrhs;

//     dmumps_c(&mumps);

//     if (mumps.info[0] != 0) {
//         std::cerr << "Error solving the system: " << mumps.info[0] << std::endl;
//         return 1;
//     }

//     std::cout << "Solution:" << std::endl;
//     for (int i = 0; i < 2*n; ++i) {
//         std::cout << mumps.rhs[i] << std::endl;
//     }


//     return 0;
// }


//     id.job = -1; // Initialize MUMPS
//     dmumps_c(&id);
    
//     // Check if MUMPS initialization was successful
//     if (id.ICNTL(1) < 0) {
//         printf("MUMPS initialization failed with error code %d\n", id.ICNTL(1));
//         return 1;
//     }
    
//     printf("MUMPS initialization successful!\n");

//     // Finalize MUMPS
//     id.job = -2; // Terminate MUMPS
//     dmumps_c(&id);

//     return 0;
// }

// #include <thread>

int main(int argc, char* argv[]){
    #ifdef NDEBUG
        std::cout << "Build Type: Release\n";
    #else
        std::cout << "Build Type: Debug\n";
    #endif

    // unsigned int numCores = std::thread::hardware_concurrency();
    // std::cout << "Number of cores available: " << numCores << std::endl;

    // allow refining of the grid at r_jump, the center point of the 
    // drop of the diffusion coefficient alpha.
    double r_jump = 0.5;
    alpha_coeff alpha_name = SONNENDRUCKER;
    if (alpha_name == SONNENDRUCKER) {
        // The center of the coefficient jump lies at 0.6888697651782026
        // for backward stability with previous runs and the Matlab code, 
        // we use 0.66 though.
        r_jump = 0.66;
    } else if (alpha_name == ZONI) {
        r_jump = 0.4837;
    } else if (alpha_name == ZONI_SHIFTED) {
        // Choose center point of descent.
        // a) - ln(0.5 * (alpha(0) - alpha(Rmax))):
        //    - ln(0.5 * (np.exp(-np.tanh(-14)) - np.exp(-np.tanh(6)))) = 0.16143743821247852
        // b) r_center = Rmax * (np.arctanh(0.16143743821247852) + 14) / 20 = 0.7081431124450334 Rmax
        r_jump = 0.7081;
    } else if (alpha_name == POISSON) {
        r_jump = 0.5; // There is no jump for Poisson so this is an arbitrary choice
    } else {
        throw std::runtime_error("Unknown alpha coeff");
    }

    std::shared_ptr<TransformationHelper> trafo = std::make_shared<TransformationHelper>();

    Fx_Functor Fx(trafo);
    Fy_Functor Fy(trafo);
    dFx_dr_Functor dFx_dr(trafo);
    dFy_dr_Functor dFy_dr(trafo);
    dFx_dt_Functor dFx_dt(trafo);
    dFy_dt_Functor dFy_dt(trafo);

    alpha_Functor alpha;
    beta_Functor beta;
    rhs_f_Functor rhs_f;
    u_D_Functor u_D;

    exact_solution_Functor exact_solution;

    GMGPolar solver;
    solver.setRadialRefinement(r_jump);
    solver.setGeometry(dFx_dr, dFy_dr, dFx_dt, dFy_dt);
    solver.setParameters(alpha, beta, rhs_f, u_D);
    solver.setSystemParameters(exact_solution);
    solver.setParameters(argc, argv);
    solver.setup(); 
    solver.solve(); 
}







// int main(int argc, char* argv[]){
//     #ifdef NDEBUG
//         std::cout << "Build Type: Release\n";
//     #else
//         std::cout << "Build Type: Debug\n";
//     #endif

//     GMGPolar solver;
//     solver.setParameters(argc, argv);
//     solver.setup();
//     solver.solve();

//     // ./compile Release
//     // ./build/gmgpolar --nr_exp 3 --ntheta_exp 4 --openmp 15

//     // Smoother Tasking //

//     int Circles = 300;
//     int BlackCircles = 150;
//     int WhiteCircles = 150;
//     // BlackCirclesTaskIndex [8,6,4,2,0]
//     // WhiteCircleTaskIndex [9,7,5,3,1]

//     int Radials = 200;
//     int BlackRadials = 100;
//     int WhiteRadials = 100;
//     // BlackRadialTaskIndex [10,12,14,16,18]
//     // WhiteRadialTaskIndex [11,13,15,17,19]

//     assert(BlackCircles == WhiteCircles || BlackCircles == WhiteCircles + 1);
//     std::vector<std::mutex> WhiteCircleMutexes(WhiteCircles);
//     std::vector<int> WhiteCircleDepCounter(WhiteCircles, 2);
//     if(BlackCircles == WhiteCircles) WhiteCircleDepCounter.back() = 1;


//     assert(BlackRadials == WhiteRadials);
//     std::vector<std::mutex> WhiteRadialMutexes(WhiteRadials);
//     std::vector<int> WhiteRadialDepCounter(WhiteRadials, 2);

    
//     auto start_time = std::chrono::high_resolution_clock::now();

//     #pragma omp parallel
//     #pragma omp single
//     {
//         // Black Circle Task 0
//         #pragma omp task
//         {
//             int localBlackCircleIndex = 0;
//             // std::cout<<"Start: Black Circle 0.\n";
//             // // Do some work...
//             std::this_thread::sleep_for(std::chrono::milliseconds(2));
//             // // Work finished
//             // std::cout<<"Finished: Black Circle 0.\n";

//             int localLeftWhiteCircleIndex = -1;
//             if(BlackCircles == WhiteCircles || localBlackCircleIndex < BlackCircles - 1){
//                 localLeftWhiteCircleIndex = ((localBlackCircleIndex) + WhiteCircles) % WhiteCircles;
//             }
//             if(localLeftWhiteCircleIndex != -1){
//                 WhiteCircleMutexes[localLeftWhiteCircleIndex].lock();
//                 WhiteCircleDepCounter[localLeftWhiteCircleIndex]--;
//                 bool leftWhiteCircle_resolved = WhiteCircleDepCounter[localLeftWhiteCircleIndex] == 0;
//                 WhiteCircleMutexes[localLeftWhiteCircleIndex].unlock();

//                 if(leftWhiteCircle_resolved){
//                     // White Circle Task
//                     #pragma omp task
//                     {
//                         // std::cout<<"Start: White Circle "<<localLeftWhiteCircleIndex<< ".\n";
//                         // // Do some work...
//                         std::this_thread::sleep_for(std::chrono::milliseconds(2));
//                         // // Work finished
//                         // std::cout<<"Finished: White Circle "<<localLeftWhiteCircleIndex<< ".\n";
//                     }
//                 }
//             }

//             int localRightWhiteCircleIndex = -1;
//             if(localBlackCircleIndex != 0){
//                 localRightWhiteCircleIndex = ((localBlackCircleIndex-1) + WhiteCircles) % WhiteCircles;
//             }
//             if(localRightWhiteCircleIndex != -1){
//                 WhiteCircleMutexes[localRightWhiteCircleIndex].lock();
//                 WhiteCircleDepCounter[localRightWhiteCircleIndex]--;
//                 bool rightWhiteCircle_resolved = WhiteCircleDepCounter[localRightWhiteCircleIndex] == 0;
//                 WhiteCircleMutexes[localRightWhiteCircleIndex].unlock();

//                 if(rightWhiteCircle_resolved){
//                     // White Circle Task
//                     #pragma omp task
//                     {
//                         // std::cout<<"Start: White Circle "<<localRightWhiteCircleIndex<< ".\n";
//                         // // Do some work...
//                         std::this_thread::sleep_for(std::chrono::milliseconds(2));
//                         // // Work finished
//                         // std::cout<<"Finished: White Circle "<<localRightWhiteCircleIndex<< ".\n";
//                     }
//                 }
//             }



//             for (int localBlackRadialIndex = 0; localBlackRadialIndex < BlackRadials; localBlackRadialIndex++)
//             {
//                 // Black Radial Task
//                 #pragma omp task
//                 {
//                     // std::cout<<"Start: Black Radial "<<localBlackRadialIndex<< ".\n";
//                     // // Do some work...
//                     std::this_thread::sleep_for(std::chrono::milliseconds(2));
//                     // // Work finished
//                     // std::cout<<"Finished: Black Radial "<<localBlackRadialIndex<< ".\n";

//                     int localLeftWhiteRadialIndex = ((localBlackRadialIndex-1) + WhiteRadials) % WhiteRadials;
//                     WhiteRadialMutexes[localLeftWhiteRadialIndex].lock();
//                     WhiteRadialDepCounter[localLeftWhiteRadialIndex]--;
//                     bool leftWhiteRadial_resolved = WhiteRadialDepCounter[localLeftWhiteRadialIndex] == 0;
//                     WhiteRadialMutexes[localLeftWhiteRadialIndex].unlock();

//                     if(leftWhiteRadial_resolved){
//                         // White Radial Task
//                         #pragma omp task
//                         {
//                             // std::cout<<"Start: White Radial "<<localLeftWhiteRadialIndex<< ".\n";
//                             // // Do some work...
//                             std::this_thread::sleep_for(std::chrono::milliseconds(2));
//                             // // Work finished
//                             // std::cout<<"Finished: White Radial "<<localLeftWhiteRadialIndex<< ".\n";
//                         }
//                     }

//                     int localRightWhiteRadialIndex = ((localBlackRadialIndex) + WhiteRadials) % WhiteRadials;
//                     WhiteRadialMutexes[localRightWhiteRadialIndex].lock();
//                     WhiteRadialDepCounter[localRightWhiteRadialIndex]--;
//                     bool rightWhiteRadial_resolved = WhiteRadialDepCounter[localRightWhiteRadialIndex] == 0;
//                     WhiteRadialMutexes[localRightWhiteRadialIndex].unlock();

//                     if(rightWhiteRadial_resolved){
//                         // White Radial Task
//                         #pragma omp task
//                         {
//                             // std::cout<<"Start: White Radial "<<localRightWhiteRadialIndex<< ".\n";
//                             // // Do some work...
//                             std::this_thread::sleep_for(std::chrono::milliseconds(2));
//                             // // Work finished
//                             // std::cout<<"Finished: White Radial "<<localRightWhiteRadialIndex<< ".\n";
//                         }
//                     }

//                 }
//             }
//         } // END OF BLACK CIRCLE 0

//         for (int localBlackCircleIndex = 1; localBlackCircleIndex < BlackCircles; localBlackCircleIndex++)
//         {
//             // Black Circle Task
//             #pragma omp task
//             {
//                 // std::cout<<"Start: Black Circle "<<localBlackCircleIndex<< ".\n";
//                 // // Do some work...
//                 std::this_thread::sleep_for(std::chrono::milliseconds(2));
//                 // // Work finished
//                 // std::cout<<"Finished: Black Circle "<<localBlackCircleIndex<< ".\n";

//                 int localLeftWhiteCircleIndex = -1;
//                 if(BlackCircles == WhiteCircles || localBlackCircleIndex < BlackCircles - 1){
//                     localLeftWhiteCircleIndex = ((localBlackCircleIndex) + WhiteCircles) % WhiteCircles;
//                 }
//                 if(localLeftWhiteCircleIndex != -1){
//                     WhiteCircleMutexes[localLeftWhiteCircleIndex].lock();
//                     WhiteCircleDepCounter[localLeftWhiteCircleIndex]--;
//                     bool leftWhiteCircle_resolved = WhiteCircleDepCounter[localLeftWhiteCircleIndex] == 0;
//                     WhiteCircleMutexes[localLeftWhiteCircleIndex].unlock();

//                     if(leftWhiteCircle_resolved){
//                         // White Circle Task
//                         #pragma omp task
//                         {
//                             // std::cout<<"Start: White Circle "<<localLeftWhiteCircleIndex<< ".\n";
//                             // // Do some work...
//                             std::this_thread::sleep_for(std::chrono::milliseconds(2));
//                             // // Work finished
//                             // std::cout<<"Finished: White Circle "<<localLeftWhiteCircleIndex<< ".\n";
//                         }
//                     }
//                 }

//                 int localRightWhiteCircleIndex = -1;
//                 if(localBlackCircleIndex != 0){
//                     localRightWhiteCircleIndex = ((localBlackCircleIndex-1) + WhiteCircles) % WhiteCircles;
//                 }
//                 if(localRightWhiteCircleIndex != -1){
//                     WhiteCircleMutexes[localRightWhiteCircleIndex].lock();
//                     WhiteCircleDepCounter[localRightWhiteCircleIndex]--;
//                     bool rightWhiteCircle_resolved = WhiteCircleDepCounter[localRightWhiteCircleIndex] == 0;
//                     WhiteCircleMutexes[localRightWhiteCircleIndex].unlock();

//                     if(rightWhiteCircle_resolved){
//                         // White Circle Task
//                         #pragma omp task
//                         {
//                             // std::cout<<"Start: White Circle "<<localRightWhiteCircleIndex<< ".\n";
//                             // // Do some work...
//                             std::this_thread::sleep_for(std::chrono::milliseconds(2));
//                             // // Work finished
//                             // std::cout<<"Finished: White Circle "<<localRightWhiteCircleIndex<< ".\n";
//                         }
//                     }
//                 }
//             }
//         }
//     }
    
//     std::cout<<"\nWhiteCircleDepCounter: ";
//     for (int i = 0; i < WhiteCircleDepCounter.size(); i++)
//     {
//         std::cout<<WhiteCircleDepCounter[i]<<", ";
//     }
    
//     std::cout<<"\nWhiteRadialDepCounter: ";
//     for (int i = 0; i < WhiteRadialDepCounter.size(); i++)
//     {
//         std::cout<<WhiteRadialDepCounter[i]<<", ";
//     }
    

//     auto end_time = std::chrono::high_resolution_clock::now();

//     auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
//     std::cout << "Execution time: " << duration_ms.count() << " milliseconds" << std::endl;

//     return 0;
// }





// int main(int argc, char* argv[]){
//     #ifdef NDEBUG
//         std::cout << "Build Type: Release\n";
//     #else
//         std::cout << "Build Type: Debug\n";
//     #endif

//     GMGPolar solver;
//     solver.setParameters(argc, argv);
//     solver.setup();
//     solver.solve();

//     // ./compile Release
//     // ./build/gmgpolar --nr_exp 3 --ntheta_exp 4 --openmp 15

//     // Smoother Tasking //

//     int Circles = 11;
//     int BlackCircles = 6;
//     int WhiteCircles = 5;
//     // BlackCirclesTaskIndex [8,6,4,2,0]
//     // WhiteCircleTaskIndex [9,7,5,3,1]

//     int Radials = 20;
//     int BlackRadials = 10;
//     int WhiteRadials = 10;
//     // BlackRadialTaskIndex [10,12,14,16,18]
//     // WhiteRadialTaskIndex [11,13,15,17,19]

//     assert(BlackCircles == WhiteCircles || BlackCircles == WhiteCircles + 1);
//     std::vector<std::mutex> WhiteCircleMutexes(WhiteCircles);
//     std::vector<int> WhiteCircleDepCounter(WhiteCircles, 2);
//     if(BlackCircles == WhiteCircles) WhiteCircleDepCounter.back() = 1;


//     assert(BlackRadials == WhiteRadials);
//     std::vector<std::mutex> WhiteRadialMutexes(WhiteRadials);
//     std::vector<int> WhiteRadialDepCounter(WhiteRadials, 2);

    
//     auto start_time = std::chrono::high_resolution_clock::now();

//     #pragma omp parallel
//     #pragma omp single
//     {
//         // Black Circle Task 0
//         #pragma omp task
//         {
//             int localBlackCircleIndex = 0;
//             std::cout<<"Start: Black Circle 0.\n";
//             // Do some work...
//             std::this_thread::sleep_for(std::chrono::milliseconds(1000));
//             // Work finished
//             std::cout<<"Finished: Black Circle 0.\n";

//             int localLeftWhiteCircleIndex = -1;
//             if(BlackCircles == WhiteCircles || localBlackCircleIndex < BlackCircles - 1){
//                 localLeftWhiteCircleIndex = ((localBlackCircleIndex) + WhiteCircles) % WhiteCircles;
//             }
//             if(localLeftWhiteCircleIndex != -1){
//                 WhiteCircleMutexes[localLeftWhiteCircleIndex].lock();
//                 WhiteCircleDepCounter[localLeftWhiteCircleIndex]--;
//                 bool leftWhiteCircle_resolved = WhiteCircleDepCounter[localLeftWhiteCircleIndex] == 0;
//                 WhiteCircleMutexes[localLeftWhiteCircleIndex].unlock();

//                 if(leftWhiteCircle_resolved){
//                     // White Circle Task
//                     #pragma omp task
//                     {
//                         std::cout<<"Start: White Circle "<<localLeftWhiteCircleIndex<< ".\n";
//                         // Do some work...
//                         std::this_thread::sleep_for(std::chrono::milliseconds(1000));
//                         // Work finished
//                         std::cout<<"Finished: White Circle "<<localLeftWhiteCircleIndex<< ".\n";
//                     }
//                 }
//             }

//             int localRightWhiteCircleIndex = -1;
//             if(localBlackCircleIndex != 0){
//                 localRightWhiteCircleIndex = ((localBlackCircleIndex-1) + WhiteCircles) % WhiteCircles;
//             }
//             if(localRightWhiteCircleIndex != -1){
//                 WhiteCircleMutexes[localRightWhiteCircleIndex].lock();
//                 WhiteCircleDepCounter[localRightWhiteCircleIndex]--;
//                 bool rightWhiteCircle_resolved = WhiteCircleDepCounter[localRightWhiteCircleIndex] == 0;
//                 WhiteCircleMutexes[localRightWhiteCircleIndex].unlock();

//                 if(rightWhiteCircle_resolved){
//                     // White Circle Task
//                     #pragma omp task
//                     {
//                         std::cout<<"Start: White Circle "<<localRightWhiteCircleIndex<< ".\n";
//                         // Do some work...
//                         std::this_thread::sleep_for(std::chrono::milliseconds(1000));
//                         // Work finished
//                         std::cout<<"Finished: White Circle "<<localRightWhiteCircleIndex<< ".\n";
//                     }
//                 }
//             }



//             for (int localBlackRadialIndex = 0; localBlackRadialIndex < BlackRadials; localBlackRadialIndex++)
//             {
//                 // Black Radial Task
//                 #pragma omp task
//                 {
//                     std::cout<<"Start: Black Radial "<<localBlackRadialIndex<< ".\n";
//                     // Do some work...
//                     std::this_thread::sleep_for(std::chrono::milliseconds(500));
//                     // Work finished
//                     std::cout<<"Finished: Black Radial "<<localBlackRadialIndex<< ".\n";

//                     int localLeftWhiteRadialIndex = ((localBlackRadialIndex-1) + WhiteRadials) % WhiteRadials;
//                     WhiteRadialMutexes[localLeftWhiteRadialIndex].lock();
//                     WhiteRadialDepCounter[localLeftWhiteRadialIndex]--;
//                     bool leftWhiteRadial_resolved = WhiteRadialDepCounter[localLeftWhiteRadialIndex] == 0;
//                     WhiteRadialMutexes[localLeftWhiteRadialIndex].unlock();

//                     if(leftWhiteRadial_resolved){
//                         // White Radial Task
//                         #pragma omp task
//                         {
//                             std::cout<<"Start: White Radial "<<localLeftWhiteRadialIndex<< ".\n";
//                             // Do some work...
//                             std::this_thread::sleep_for(std::chrono::milliseconds(500));
//                             // Work finished
//                             std::cout<<"Finished: White Radial "<<localLeftWhiteRadialIndex<< ".\n";
//                         }
//                     }

//                     int localRightWhiteRadialIndex = ((localBlackRadialIndex) + WhiteRadials) % WhiteRadials;
//                     WhiteRadialMutexes[localRightWhiteRadialIndex].lock();
//                     WhiteRadialDepCounter[localRightWhiteRadialIndex]--;
//                     bool rightWhiteRadial_resolved = WhiteRadialDepCounter[localRightWhiteRadialIndex] == 0;
//                     WhiteRadialMutexes[localRightWhiteRadialIndex].unlock();

//                     if(rightWhiteRadial_resolved){
//                         // White Radial Task
//                         #pragma omp task
//                         {
//                             std::cout<<"Start: White Radial "<<localRightWhiteRadialIndex<< ".\n";
//                             // Do some work...
//                             std::this_thread::sleep_for(std::chrono::milliseconds(500));
//                             // Work finished
//                             std::cout<<"Finished: White Radial "<<localRightWhiteRadialIndex<< ".\n";
//                         }
//                     }

//                 }
//             }
//         } // END OF BLACK CIRCLE 0

//         for (int localBlackCircleIndex = 1; localBlackCircleIndex < BlackCircles; localBlackCircleIndex++)
//         {
//             // Black Circle Task
//             #pragma omp task
//             {
//                 std::cout<<"Start: Black Circle "<<localBlackCircleIndex<< ".\n";
//                 // Do some work...
//                 std::this_thread::sleep_for(std::chrono::milliseconds(1000));
//                 // Work finished
//                 std::cout<<"Finished: Black Circle "<<localBlackCircleIndex<< ".\n";

//                 int localLeftWhiteCircleIndex = -1;
//                 if(BlackCircles == WhiteCircles || localBlackCircleIndex < BlackCircles - 1){
//                     localLeftWhiteCircleIndex = ((localBlackCircleIndex) + WhiteCircles) % WhiteCircles;
//                 }
//                 if(localLeftWhiteCircleIndex != -1){
//                     WhiteCircleMutexes[localLeftWhiteCircleIndex].lock();
//                     WhiteCircleDepCounter[localLeftWhiteCircleIndex]--;
//                     bool leftWhiteCircle_resolved = WhiteCircleDepCounter[localLeftWhiteCircleIndex] == 0;
//                     WhiteCircleMutexes[localLeftWhiteCircleIndex].unlock();

//                     if(leftWhiteCircle_resolved){
//                         // White Circle Task
//                         #pragma omp task
//                         {
//                             std::cout<<"Start: White Circle "<<localLeftWhiteCircleIndex<< ".\n";
//                             // Do some work...
//                             std::this_thread::sleep_for(std::chrono::milliseconds(1000));
//                             // Work finished
//                             std::cout<<"Finished: White Circle "<<localLeftWhiteCircleIndex<< ".\n";
//                         }
//                     }
//                 }

//                 int localRightWhiteCircleIndex = -1;
//                 if(localBlackCircleIndex != 0){
//                     localRightWhiteCircleIndex = ((localBlackCircleIndex-1) + WhiteCircles) % WhiteCircles;
//                 }
//                 if(localRightWhiteCircleIndex != -1){
//                     WhiteCircleMutexes[localRightWhiteCircleIndex].lock();
//                     WhiteCircleDepCounter[localRightWhiteCircleIndex]--;
//                     bool rightWhiteCircle_resolved = WhiteCircleDepCounter[localRightWhiteCircleIndex] == 0;
//                     WhiteCircleMutexes[localRightWhiteCircleIndex].unlock();

//                     if(rightWhiteCircle_resolved){
//                         // White Circle Task
//                         #pragma omp task
//                         {
//                             std::cout<<"Start: White Circle "<<localRightWhiteCircleIndex<< ".\n";
//                             // Do some work...
//                             std::this_thread::sleep_for(std::chrono::milliseconds(1000));
//                             // Work finished
//                             std::cout<<"Finished: White Circle "<<localRightWhiteCircleIndex<< ".\n";
//                         }
//                     }
//                 }
//             }
//         }
//     }
    
//     std::cout<<"\nWhiteCircleDepCounter: ";
//     for (int i = 0; i < WhiteCircleDepCounter.size(); i++)
//     {
//         std::cout<<WhiteCircleDepCounter[i]<<", ";
//     }
    
//     std::cout<<"\nWhiteRadialDepCounter: ";
//     for (int i = 0; i < WhiteRadialDepCounter.size(); i++)
//     {
//         std::cout<<WhiteRadialDepCounter[i]<<", ";
//     }
    

//     auto end_time = std::chrono::high_resolution_clock::now();

//     auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
//     std::cout << "Execution time: " << duration_ms.count() << " milliseconds" << std::endl;

//     return 0;
// }






        //     for (int localBlackCircleIndex = 0; localBlackCircleIndex < WhiteCircles; localBlackCircleIndex++)
        //     {
        //         // BlackCircle Task localBlackCircleIndex
        //         #pragma omp task
        //         {
        //             std::cout<<"Start: Black Circle "<<localBlackCircleIndex<< ".\n";
        //             // Do some work...
        //             std::this_thread::sleep_for(std::chrono::milliseconds(10000));
        //             // Work finished
        //             std::cout<<"Finished: Black Circle "<<localBlackCircleIndex<< ".\n";

        //             WhiteCircleMutexes[0].lock();

        //             WhiteCircleDepCounter[0]--;
        //             bool resolved = WhiteCircleDepCounter[0] == 0;

        //             WhiteCircleMutexes[0].unlock();

        //             if(resolved){

        //                 // WhiteCircle Task localWhiteCircleIndex
        //                 int localWhiteCircleIndex = 
        //                 #pragma omp task
        //                 {
        //                     std::cout<<"Start: Black Circle "<<localWhiteCircleIndex<< ".\n";
        //                     // Do some work...
        //                     std::this_thread::sleep_for(std::chrono::milliseconds(10000));
        //                     // Work finished
        //                     std::cout<<"Finished: Black Circle "<<localWhiteCircleIndex<< ".\n";

        //                     }


        //         }
        //     }
        // }




// chmod +x compile.sh
// ./compile.sh Release
// ./build/gmgpolar --openmp 10
// ./build/gmgpolar --openmp 50
// ./build/gmgpolar --openmp 100

// void executeTasks() {
//     // Tasks. Applys also to small values of T.
//     const int T = 30'000;
//     int* dep = new int[T];

//     double start = omp_get_wtime();

//     omp_set_dynamic(1);

//     #pragma omp parallel shared(dep)
//     {
//         #pragma omp single
//         {
//             for(int i = 0; i < T; i += 3) {
//                 #pragma omp task shared(dep) firstprivate(i) depend(out: dep[i])
//                 {
//                 }
//             }
//         }
//     }

//     double end = omp_get_wtime(); // End time
//     double duration = (end - start) * 1000.0; // Calculate duration in milliseconds

//     std::cout << "Time taken for tasks: " << duration << " milliseconds" << std::endl;

//     delete[] dep;
// }

// int main(int argc, char* argv[]){
//     #ifdef NDEBUG
//         std::cout << "Build Type: Release\n";
//     #else
//         std::cout << "Build Type: Debug\n";
//     #endif

//     const int numThreads = omp_get_max_threads();
//     std::cout<<"Threads available: "<<numThreads<<std::endl;

//     GMGPolar solver;
//     solver.setParameters(argc, argv); // Just here for setting openMP threads

//     const int usedThreads = omp_get_max_threads();
//     std::cout<<"Threads used: "<<usedThreads<<std::endl;

//     executeTasks();
//     return 0;
// }



// int main(int argc, char* argv[]){
//     #ifdef NDEBUG
//         std::cout << "Build Type: Release\n";
//     #else
//         std::cout << "Build Type: Debug\n";
//     #endif

//     // int numTasks = 2;
//     // int minChunkSize = 3;
//     // int numThreads = 4;
//     // TaskDistribution(numTasks,minChunkSize,numThreads);

//     GMGPolar solver;
//     solver.setParameters(argc, argv);
//     solver.setup();
//     return 0;
// }

// #include <iostream>
// #include <vector>
// #include <omp.h>
// #include <thread>
// #include <chrono>
// #include <random>

// void executeTasks() {
//     const int numThreads = omp_get_max_threads();
//     std::cout<<"Used: "<<numThreads<<std::endl;
//     omp_set_num_threads(50);

//     const int S1 = 100;
//     const int S2 = 67;
//     const int T = S1 + S2;

//     int* dep = new int[T];

//     const int S2_wait = S2 % 3;
//     const int S2_start = std::max(S1-2, 0);

//     #pragma omp parallel shared(dep)
//     {
//         #pragma omp single
//         {
//             for(int i = S1 - 1; i >= 0; i -= 3) {
//                 #pragma omp task shared(dep) firstprivate(i) depend(out: dep[i])
//                 {
//                     std::cout << "Task " << i << " finished" << std::endl;
//                 }
//             }

//             for(int i = S1 - 2; i >= 0; i -= 3) {
//                 #pragma omp task shared(dep) firstprivate(i) depend(out: dep[i]) depend(in: dep[i+1], dep[i-2])
//                 {
//                     std::cout << "Task " << i << " finished" << std::endl;
//                 }
//             }

//             for(int i = S1 - 3; i >= 0; i -= 3) {
//                 #pragma omp task shared(dep) firstprivate(i) depend(out: dep[i]) depend(in: dep[i+1], dep[i-2])
//                 {
//                     std::cout << "Task " << i << " finished" << std::endl;
//                 }
//             }

//             for(int i = S1; i < T - S2_wait; i += 3) {
//                 #pragma omp task shared(dep) firstprivate(i) depend(out: dep[i]) depend(in: dep[S2_start])
//                 {
//                     std::cout << "Task " << i << " finished" << std::endl;
//                 }
//             }

//             for(int i = S1+1; i < T - S2_wait; i += 3) {
//                 #pragma omp task shared(dep) firstprivate(i) depend(out: dep[i]) depend(in: dep[S2_start],  dep[i-1], dep[i+2])
//                 {
//                     std::cout << "Task " << i << " finished" << std::endl;
//                 }
//             }

//             for(int i = S1+2; i < T - S2_wait; i += 3) {
//                 #pragma omp task shared(dep) firstprivate(i) depend(out: dep[i]) depend(in: dep[S2_start],  dep[i-1], dep[i+2])
//                 {
//                     std::cout << "Task " << i << " finished" << std::endl;
//                 }
//             }

//             if( S2_wait >= 1 ){
//                 int i = T - S2_wait;
//                 #pragma omp task shared(dep) firstprivate(i) depend(out: dep[i]) depend(in: dep[S2_start],  dep[i-1])
//                 {
//                     std::cout << "Task " << i << " finished" << std::endl;
//                 }
//             }

//             if( S2_wait >= 2 ){
//                 int i = T - S2_wait + 1;
//                 #pragma omp task shared(dep) firstprivate(i) depend(out: dep[i]) depend(in: dep[S2_start],  dep[i-1])
//                 {
//                     std::cout << "Task " << i << " finished" << std::endl;
//                 }
//             }
//         }
//     }

//     delete[] dep;


//     // // Execute tasks
//     // #pragma omp parallel shared(dep)
//     // {
//     //     #pragma omp single
//     //     {
//     //         // Generate tasks with correct dependencies
//     //         for (int i : loop) {
//     //             // Create a separate depend clause for each dependency
//     //             #pragma omp task shared(dep) firstprivate(i) depend(in: dep[dependencies[k] : ])



//     //             for(int k = 0; k < dependencies[i].size(); k++) {
//     //                 #pragma omp depend(in: dep[dependencies[k]])
//     //             }
//     //             {
//     //                 // Task execution
//     //                     std::cout << "Executing Task " << i << std::endl;    
//     //                     std::this_thread::sleep_for(std::chrono::duration<double>(0.1)); // Wait for the random time

//     //                 // Mark task as completed
//     //                 dep[i] = 1;
//     //             }
//     //         }
//     //     }
//     // }
// }

// int main() {
//     executeTasks();
//     return 0;
// }

















// void level::build_A()
// {
//     int start_j;
//     int* dep = new int[nr];

//     // Take boundary condition into account: Dirichlet-RB
//     if (gyro::icntl[Param::DirBC_Interior]) { // (r[0],0) is on Dirichlet boundary
//         start_j = 1;
//     }
//     else { // (r[0],0) is not on Dirichlet boundary
//         start_j = 0;
//     }

//     #pragma omp parallel shared(dep)
//     {
//         #pragma omp single
//         {
//             if (gyro::icntl[Param::DirBC_Interior]) { // (r[0],0) is on Dirichlet boundary
//                 #pragma omp task shared(dep, start_j) depend(out : dep[0])
//                 {
//                     for (int i = 0; i < ntheta_int; i++) {}
//                 } //end of task and parallel
//             }
//             #pragma omp task shared(dep, start_j) depend(out : dep[nr_int])
//             {
//                 // Take boundary condition into account: Dirichlet-RB
//                 for (int i = 0; i < ntheta_int; i++){}
//             } //end of task and parallel

//             #pragma omp task shared(dep, start_j) depend(out : dep[start_j])
//             {
//                 // ------------ //
//                 // Inner Circle //
//                 // ------------ //

//                 // center
//                 for (int i = 0; i < ntheta_int; i++) {
//                     vals[ptr + stencil[Param::middle]] += val;
//                 }

//                 // Across: bottom update
//                 if (!gyro::icntl[Param::DirBC_Interior]) {
//                     for (i = 0; i < ntheta_int; i++) {
//                         vals[ptr + stencil[Param::bottom]] += -coeff - coeff2;
//                         vals[ptr + stencil[Param::middle]] += coeff;
//                     }
//                 }
//                 else {
//                     for (i = 0; i < ntheta_int; i++) {
//                         vals[ptr + stencil[Param::middle]] += coeff;
//                     }
//                 }

//                 // Across and DB_int updates (~~~ Interior - (j-1, i))
//                 for (i = 0; i < ntheta_int; i++) {
//                     vals[ptr + stencil[Param::right]] += -coeff2 / kt;
//                     vals[ptr + stencil[Param::middle]] += coeff / hs + coeff / hsmin1 + coeff2 / ktmin1 + coeff2 / kt;
//                     coeff = 0.5 * (hs + hsmin1) * att_vect[i] / kt;
//                     vals[ptr + stencil[Param::left]] += -coeff;
//                     vals[ptr + stencil[Param::middle]] += coeff;

//                     // Update (j, i-1) (=== Interior - bottom_right)
//                     coeff = 0.5 * (hs + hsmin1) * att_vect[i] / ktmin1;
//                     vals[ptr + stencil[Param::right]] += -coeff;
//                     vals[ptr + stencil[Param::middle]] += coeff;

//                     // Update (j+1, i) (=== Interior)
//                     coeff = 0.5 * (kt + ktmin1) * arr_vect[i] / hs;
//                     vals[ptr + stencil[Param::bottom]] += -coeff;
//                     vals[ptr + stencil[Param::middle]] += coeff;
//                 }
//                 if (gyro::icntl[Param::mod_pk] > 0)
//                     for (i = 0; i < ntheta_int; i++) {
//                         vals[ptr + stencil[Param::top_left]] += art_vect[i];
//                         vals[ptr + stencil[Param::top_right]] += -art_vect[i];
//                         // Update (j+1, i) (=== Interior)
//                         vals[ptr + stencil[Param::bottom_left]] += -art_vect[i];
//                         vals[ptr + stencil[Param::bottom_right]] += art_vect[i];
//                     }
//             } // end of task

//             // ---------- //
//             // Mod 0 Line //
//             // ---------- //
//             for (int j = start_j + 3; j < nr_int - 1; j += 3) {
//                 #pragma omp task shared(dep, start_j) firstprivate(j) depend(out : dep[j])
//                 {
//                     for (int i = 0; i < ntheta_int; i++) {
//                         vals[ptr + stencil[Param::middle]] += val;
//                     }

//                     for (int i = 0; i < ntheta_int; i++) {
//                         coeff = 0.5 * (kt + ktmin1) * arr_vect[i] / hsmin1;
//                         vals[ptr + stencil[Param::top]] += -coeff;
//                         vals[ptr + stencil[Param::middle]] += coeff;

//                         // Update (j, i)
//                         coeff = 0.5 * (kt + ktmin1) * arr_vect[i];
//                         vals[ptr + stencil[Param::top]] += -coeff / hs;
//                         vals[ptr + stencil[Param::bottom]] += -coeff / hsmin1;

//                         coeff2 = 0.5 * (hs + hsmin1) * att_vect[i];
//                         vals[ptr + stencil[Param::left]] += -coeff2 / ktmin1;
//                         vals[ptr + stencil[Param::right]] += -coeff2 / kt;
//                         vals[ptr + stencil[Param::middle]] +=
//                             coeff / hs + coeff / hsmin1 + coeff2 / ktmin1 + coeff2 / kt;

//                         // Update (j, i+1)
//                         coeff = 0.5 * (hs + hsmin1) * att_vect[i] / kt;
//                         vals[ptr + stencil[Param::left]] += -coeff;
//                         vals[ptr + stencil[Param::middle]] += coeff;

//                         // Update (j, i-1)
//                         coeff = 0.5 * (hs + hsmin1) * att_vect[i] / ktmin1;
//                         vals[ptr + stencil[Param::right]] += -coeff;
//                         vals[ptr + stencil[Param::middle]] += coeff;

//                         // Update (j+1, i) (Not in DB_ext)
//                         coeff = 0.5 * (kt + ktmin1) * arr_vect[i] / hs;
//                         vals[ptr + stencil[Param::bottom]] += -coeff;
//                         vals[ptr + stencil[Param::middle]] += coeff;
//                     }
//                     if (gyro::icntl[Param::mod_pk] > 0)
//                         for (int i = 0; i < ntheta_int; i++) {
//                             // Update (j-1, i) (Not in DB_int)
//                             vals[ptr + stencil[Param::top_left]] += art_vect[i];
//                             vals[ptr + stencil[Param::top_right]] += -art_vect[i];
//                             // Update (j, i+1)
//                             vals[ptr + stencil[Param::bottom_left]] += -art_vect[i];
//                             vals[ptr + stencil[Param::top_left]] += art_vect[i];
//                             // Update (j, i-1)
//                             vals[ptr + stencil[Param::bottom_right]] += art_vect[i];
//                             vals[ptr + stencil[Param::top_right]] += -art_vect[i];
//                             // Update (j+1, i) (Not in DB_ext)
//                             vals[ptr + stencil[Param::bottom_left]] += -art_vect[i];
//                             vals[ptr + stencil[Param::bottom_right]] += art_vect[i];
//                         }
//                 } // end of task
//             }
//             // ---------- //
//             // Mod 1 Line //
//             // ---------- //
//             for (int j = start_j + 1; j < nr_int - 1; j += 3) {
//                 #pragma omp task shared(dep, start_j) firstprivate(j) depend(in : dep[j - 1]) depend(in : dep[j + 2]) depend(out : dep[j])
//                 {
//                     for (int i = 0; i < ntheta_int; i++) {
//                         vals[ptr + stencil[Param::middle]] += val;
//                     }
//                     for (int i = 0; i < ntheta_int; i++) {
//                         coeff = 0.5 * (kt + ktmin1) * arr_vect[i] / hsmin1;
//                         vals[ptr + stencil[Param::top]] += -coeff;
//                         vals[ptr + stencil[Param::middle]] += coeff;
//                         // Update (j, i)
//                         coeff = 0.5 * (kt + ktmin1) * arr_vect[i];
//                         vals[ptr + stencil[Param::top]] += -coeff / hs;
//                         vals[ptr + stencil[Param::bottom]] += -coeff / hsmin1;

//                         coeff2 = 0.5 * (hs + hsmin1) * att_vect[i];
//                         vals[ptr + stencil[Param::left]] += -coeff2 / ktmin1;
//                         vals[ptr + stencil[Param::right]] += -coeff2 / kt;
//                         vals[ptr + stencil[Param::middle]] +=
//                             coeff / hs + coeff / hsmin1 + coeff2 / ktmin1 + coeff2 / kt;
//                         // Update (j, i+1)
//                         coeff = 0.5 * (hs + hsmin1) * att_vect[i] / kt;
//                         vals[ptr + stencil[Param::left]] += -coeff;
//                         vals[ptr + stencil[Param::middle]] += coeff;
//                         // Update (j, i-1)
//                         coeff = 0.5 * (hs + hsmin1) * att_vect[i] / ktmin1;
//                         vals[ptr + stencil[Param::right]] += -coeff;
//                         vals[ptr + stencil[Param::middle]] += coeff;
//                         // Update (j+1, i) (Not in DB_ext)
//                         coeff = 0.5 * (kt + ktmin1) * arr_vect[i] / hs;
//                         vals[ptr + stencil[Param::bottom]] += -coeff;
//                         vals[ptr + stencil[Param::middle]] += coeff;
//                     }
//                     if (gyro::icntl[Param::mod_pk] > 0)
//                         for (int i = 0; i < ntheta_int; i++) {
//                             // Update (j-1, i) (Not in DB_int)
//                             row_indices[ptr + stencil[Param::top_left]] = row;
//                             col_indices[ptr + stencil[Param::top_left]] = col;

//                             vals[ptr + stencil[Param::top_left]] += art_vect[i];
//                             col = j * ntheta_int + ((i + 1 > ntheta_int - 1) ? 0 : i + 1);
//                             row_indices[ptr + stencil[Param::top_right]] = row;
//                             col_indices[ptr + stencil[Param::top_right]] = col;

//                             vals[ptr + stencil[Param::top_right]] += -art_vect[i];

//                             // Update (j, i+1)
//                             ptr     = ptr_vect[i + 2];
//                             stencil = stencil_cur;
//                             row     = j * ntheta_int + ((i + 1 > ntheta_int - 1) ? 0 : i + 1);
//                             col     = (j - 1) * ntheta_int + i;
//                             row_indices[ptr + stencil[Param::bottom_left]] = row;
//                             col_indices[ptr + stencil[Param::bottom_left]] = col;

//                             vals[ptr + stencil[Param::bottom_left]] += -art_vect[i];
//                             col                                         = (j + 1) * ntheta_int + i;
//                             row_indices[ptr + stencil[Param::top_left]] = row;
//                             col_indices[ptr + stencil[Param::top_left]] = col;

//                             vals[ptr + stencil[Param::top_left]] += art_vect[i];

//                             // Update (j, i-1)
//                             ptr = ptr_vect[i];
//                             row = j * ntheta_int + ((i - 1 < 0) ? ntheta_int - 1 : i - 1);
//                             col = (j - 1) * ntheta_int + i;
//                             row_indices[ptr + stencil[Param::bottom_right]] = row;
//                             col_indices[ptr + stencil[Param::bottom_right]] = col;

//                             vals[ptr + stencil[Param::bottom_right]] += art_vect[i];
//                             col                                          = (j + 1) * ntheta_int + i;
//                             row_indices[ptr + stencil[Param::top_right]] = row;
//                             col_indices[ptr + stencil[Param::top_right]] = col;

//                             vals[ptr + stencil[Param::top_right]] += -art_vect[i];

//                             // Update (j+1, i) (Not in DB_ext)
//                             ptr     = ptr_vect_next[i + 1];
//                             stencil = stencil_next;
//                             row     = (j + 1) * ntheta_int + i;
//                             col     = j * ntheta_int + ((i - 1 < 0) ? ntheta_int - 1 : i - 1);
//                             row_indices[ptr + stencil[Param::bottom_left]] = row;
//                             col_indices[ptr + stencil[Param::bottom_left]] = col;

//                             vals[ptr + stencil[Param::bottom_left]] += -art_vect[i];
//                             col = j * ntheta_int + ((i + 1 > ntheta_int - 1) ? 0 : i + 1);
//                             row_indices[ptr + stencil[Param::bottom_right]] = row;
//                             col_indices[ptr + stencil[Param::bottom_right]] = col;

//                             vals[ptr + stencil[Param::bottom_right]] += art_vect[i];
//                         }
//                 } // end of task
//             }
//             // ---------- //
//             // Mod 2 Line //
//             // ---------- //
//             for (int j = start_j + 2; j < nr_int - 1; j += 3) {
//                 #pragma omp task shared(dep, start_j) firstprivate(j) depend(in : dep[j - 1]) depend(in : dep[j + 2]) depend(out : dep[j])
//                 {
//                     for (int i = 0; i < ntheta_int; i++) {
//                         vals[ptr + stencil[Param::middle]] += val;
//                     }

//                     for (int i = 0; i < ntheta_int; i++) {
//                         // Update (j-1, i) (Not in DB_int)
//                         coeff = 0.5 * (kt + ktmin1) * arr_vect[i] / hsmin1;
//                         vals[ptr + stencil[Param::top]] += -coeff;
//                         vals[ptr + stencil[Param::middle]] += coeff;
//                         // Update (j, i)
//                         coeff = 0.5 * (kt + ktmin1) * arr_vect[i];
//                         vals[ptr + stencil[Param::top]] += -coeff / hs;
//                         vals[ptr + stencil[Param::bottom]] += -coeff / hsmin1;
//                         coeff2 = 0.5 * (hs + hsmin1) * att_vect[i];
//                         vals[ptr + stencil[Param::left]] += -coeff2 / ktmin1;
//                         vals[ptr + stencil[Param::right]] += -coeff2 / kt;
//                         vals[ptr + stencil[Param::middle]] +=
//                             coeff / hs + coeff / hsmin1 + coeff2 / ktmin1 + coeff2 / kt;
//                         // Update (j, i+1)
//                         coeff = 0.5 * (hs + hsmin1) * att_vect[i] / kt;
//                         vals[ptr + stencil[Param::left]] += -coeff;
//                         vals[ptr + stencil[Param::middle]] += coeff;
//                         // Update (j, i-1)
//                         coeff = 0.5 * (hs + hsmin1) * att_vect[i] / ktmin1;
//                         vals[ptr + stencil[Param::right]] += -coeff;
//                         vals[ptr + stencil[Param::middle]] += coeff;
//                         // Update (j+1, i) (Not in DB_ext)
//                         coeff = 0.5 * (kt + ktmin1) * arr_vect[i] / hs;
//                         vals[ptr + stencil[Param::bottom]] += -coeff;
//                         vals[ptr + stencil[Param::middle]] += coeff;
//                     }
//                     if (gyro::icntl[Param::mod_pk] > 0)
//                         for (int i = 0; i < ntheta_int; i++) {
//                             // Update (j-1, i) (Not in DB_int)
//                             vals[ptr + stencil[Param::top_left]] += art_vect[i];
//                             vals[ptr + stencil[Param::top_right]] += -art_vect[i];
//                             // Update (j, i+1)
//                             vals[ptr + stencil[Param::bottom_left]] += -art_vect[i];
//                             vals[ptr + stencil[Param::top_left]] += art_vect[i];
//                             // Update (j, i-1)
//                             vals[ptr + stencil[Param::bottom_right]] += art_vect[i];
//                             vals[ptr + stencil[Param::top_right]] += -art_vect[i];
//                             // Update (j+1, i) (Not in DB_ext)
//                             vals[ptr + stencil[Param::bottom_left]] += -art_vect[i];
//                             vals[ptr + stencil[Param::bottom_right]] += art_vect[i];
//                         }
//                 } // end of task
//             }

//             #pragma omp task shared(dep, start_j) depend(in : dep[nr_int - 2]) depend(in : dep[nr_int - 3]) depend(out : dep[nr_int - 1])
//             {
//             // ------------ //
//             // Outer Circle //
//             // ------------ //
//                 // DB_ext updates (~~~ Interior - (j+1, i) + DB)
//                 for (int i = 0; i < ntheta_int; i++) {
//                     vals[ptr + stencil[Param::middle]] += val;
//                 }

//                 for (int i = 0; i < ntheta_int; i++) {
//                     // Update (j-1, i) (=== Interior)
//                     coeff = 0.5 * (kt + ktmin1) * arr_vect[i] / hsmin1;
//                     vals[ptr + stencil[Param::top]] += -coeff;
//                     vals[ptr + stencil[Param::middle]] += coeff;
//                     // Update (j, i) (=== Interior - top)
//                     vals[ptr + stencil[Param::bottom]] += -coeff / hsmin1;
//                     // Contribution to middle (top) from DB
//                     coeff2 = 0.5 * (kt + ktmin1) * arr_vect2[i] / hs;
//                     coeff3 = 0.5 * (hs + hsmin1) * att_vect[i];
//                     vals[ptr + stencil[Param::left]] += -coeff3 / ktmin1;
//                     vals[ptr + stencil[Param::right]] += -coeff3 / kt;
//                     vals[ptr + stencil[Param::middle]] +=
//                         coeff / hsmin1 + coeff / hs + coeff2 + coeff3 / ktmin1 + coeff3 / kt;
//                     // Update (j, i+1) (=== Interior - top_left)
//                     coeff = 0.5 * (hs + hsmin1) * att_vect[i] / kt;
//                     vals[ptr + stencil[Param::left]] += -coeff;
//                     vals[ptr + stencil[Param::middle]] += coeff;
//                     // Update (j, i-1) (=== Interior - top_right)
//                     coeff = 0.5 * (hs + hsmin1) * att_vect[i] / ktmin1;
//                     vals[ptr + stencil[Param::right]] += -coeff;
//                     vals[ptr + stencil[Param::middle]] += coeff;
//                 }
//                 if (gyro::icntl[Param::mod_pk] > 0)
//                     for (int i = 0; i < ntheta_int; i++) {
//                         // Update (j-1, i) (=== Interior)
//                         vals[ptr + stencil[Param::top_left]] += art_vect[i];
//                         vals[ptr + stencil[Param::top_right]] += -art_vect[i];
//                         // Update (j, i+1) (=== Interior - top_left)
//                         vals[ptr + stencil[Param::bottom_left]] += -art_vect[i];
//                         // Update (j, i-1) (=== Interior - top_right)
//                         vals[ptr + stencil[Param::bottom_right]] += art_vect[i];
//                     }
//             } //end of task
//         } //end of single
//     } //end of parallel

//     delete[] dep;
// } /* ----- end of level::build_A ----- */
