#include <iostream>
#include <vector>
#include <cmath>
#include <omp.h>

#include "../include/GMGPolar/gmgpolar.h"

// #include "../include/TaskDistribution/TaskDistribution.h"


// #include "../include/linear_algebra/vector.h"
// #include "../include/linear_algebra/operations.h"
// #include <chrono>

// #include <iostream>
// #include <vector>
// #include <omp.h>
// #include <thread>
// #include <chrono>
// #include <random>

int main(int argc, char* argv[]){
    #ifdef NDEBUG
        std::cout << "Build Type: Release\n";
    #else
        std::cout << "Build Type: Debug\n";
    #endif

    // int numTasks = 2;
    // int minChunkSize = 3;
    // int numThreads = 4;
    // TaskDistribution(numTasks,minChunkSize,numThreads);

    GMGPolar solver;
    solver.setParameters(argc, argv);
    solver.setup();
    solver.solve();

    // std::unique_ptr<ExactFunctions> functions = std::make_unique<CartesianR2GyroSonnendruckerCircular>();

    return 0;
}




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
